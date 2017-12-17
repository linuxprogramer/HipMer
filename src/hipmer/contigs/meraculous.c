#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <assert.h>
#include <upc_tick.h>
#include <libgen.h>
#include "../common/optlist.h"
#include "../common/common.h"

#include "meraculous.h"


#ifdef LSH

shared[1] oracle_t *oracle_table;

#endif


#ifdef LHS_PERF
int64_t lookups = 0;
int64_t success = 0;
#endif


#ifdef USE_CRAY_UPC
#include <intrinsics.h>
#include <upc.h>
#endif

#ifdef USE_BUPC
#include <upc.h>
#endif

#ifdef USE_CRAY_UPC
#include <upc_cray.h>
#endif

#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include "../contigs/kmer_hash.h"
#include "kmer_handling.h"

#ifndef STORE_OPT
#include "buildUFXhash.h"
#endif

#ifdef STORE_OPT

double fileIOTime = 0.0;
double cardCalcTime = 0.0;
double setupTime = 0.0;
double storeTime = 0.0;

#ifndef BIN_INPUT
#include "buildUFXhashOPT.h"
#endif

#ifdef BIN_INPUT
#include "../kmercount/readufx.h"
#include "buildUFXhashBinary.h"
#endif

#endif

int myContigs = 0;

#include "UU_traversal_final.h"

#ifdef USE_BUPC
shared int64_t contig_id;
#endif

#ifdef USE_CRAY_UPC
shared long contig_id;
#endif
shared[BS] contig_ptr_t *contig_list;

shared int64_t timestamp;
int64_t local_allocs = 0;
 
#ifdef PROFILE
int64_t bases_added = 0;
#endif

#ifdef UU_TRAV_PROFILE
double walking_time = 0.0;
double UU_time = 0.0;
#endif

int cores_per_node = 24;

int main(int argc, char **argv) {
    upc_tick_t start_time = upc_ticks_now();
	
   int fileNo = 0;
   int i;
	int64_t size;
	hash_table_t *dist_hashtable;
	memory_heap_t memory_heap;
	FILE *fd, *blastresFD;
   FILE *my_out_file;
   FILE *reads_fd, *fdReads;
   char readFileName[MAX_FILE_PATH];
   char outputBlastmap[MAX_FILE_PATH];
   const char *base_dir = ".";
   shared[1] double *setupTimes;
   shared[1] double *calcTimes;
   shared[1] double *ioTimes;
   shared[1] double *storeTimes;
   double min_IO, min_calc, min_store, min_setup;
   double max_IO, max_store, max_calc, max_setup;
   double sum_ios=0.0;
   double sum_calcs=0.0;
   double sum_stores=0.0;
   double sum_setups=0.0;
   
   
   /* Use getopt() to read arguments */
   extern char *optarg;
	extern int optind;
   int c;
   char *input_UFX_name, *output_name, *read_files_name;
   int minimum_contig_length=MINIMUM_CONTIG_SIZE, dmin=10;
   int chunk_size = 1;
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "i:o:m:d:c:s:l:f:F:C:B:N:");
   
   char *string_size;
   int load_factor = 1;
   char *oracle_file;
   
   /* Sizes should be given in megabytes */
   int64_t outputSizeSimple;
   int64_t outputSizeFasta;
   
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'i':
            input_UFX_name = thisOpt->argument;
            break;
         case 'o':
            output_name = thisOpt->argument;
            break;
         case 'm':
            minimum_contig_length = atoi(thisOpt->argument);
            break;
         case 'd':
            dmin = atoi(thisOpt->argument);
            break;
         case 'c':
            chunk_size = atoi(thisOpt->argument);
            break;
         case 's':
            string_size = thisOpt->argument;
            break;
         case 'l':
            load_factor = atoi(thisOpt->argument);
            break;
         case 'f':
            oracle_file = thisOpt->argument;
            break;
         case 'F':
            outputSizeFasta = atoi(thisOpt->argument);
            break;
         case 'C':
            outputSizeSimple = atoi(thisOpt->argument);
            break;
         case 'B':
             base_dir = thisOpt->argument;
			 break;
         case 'N':
			 cores_per_node = atoi(thisOpt->argument);
			 break;
         default:
            break;
      }
      
      free(thisOpt);
   }
	
	UPC_TICK_T start, end;

	if (MYTHREAD == 0) {
#ifndef NO_PAD
		printf("Struct size is %lu, with %lu padding, %lu shared[] ptr\n", 
               (unsigned long)sizeof(list_t), (unsigned long)UPC_PADDING, (unsigned long)sizeof(shared void*));
#else
		printf("Struct size is %lu, no padding, %lu shared[] ptr\n", 
               (unsigned long)sizeof(list_t), (unsigned long)sizeof(shared void*));
#endif
		printf("Compiled with KMER_LENGTH %d\n", KMER_LENGTH);
		printf("Running with %d cores per node\n", cores_per_node);
	}
	if (strlen(output_name) > 200) {
		if (MYTHREAD == 0) printf("Please choose a shorter output name max 200 chars\n");
		upc_global_exit(1);
	}
	double con_time, trav_time, blastmap_time;

#if VERBOSE > 0
	{
		char fname[MAX_FILE_PATH];
		sprintf(fname, "meraculous_%s_%d.log", output_name, MYTHREAD);	
		mylog = fopen_rank_path(fname, "w", MYTHREAD);
	}
#else
	mylog = NULL;
#endif

#ifdef SYNC_PROTOCOL
	if (MYTHREAD == 0)
		printf("\n\n*************** RUNNING VERSION WITH SYNCHRONIZATION PROTOCOL ENABLED ***************\n\n");
#endif
	
#ifndef SYNC_PROTOCOL
	if (MYTHREAD == 0)
		printf("\n\n*************** RUNNING VERSION WITH SYNCHRONIZATION PROTOCOL DISABLED ***************\n\n");
#endif

#ifdef DEBUG
	if (MYTHREAD == 0)
		printf("\n\n*************** RUNNING VERSION WITH DEBUG assertions and extra checks ***************\n\n");
#endif

#if VERBOSE > 0
	if (MYTHREAD == 0)
		printf("\n\n*************** RUNNING LOGGED VERBOSE (%d) OUTPUT ***********************************\n\n", VERBOSE);
#endif


#ifdef PROFILE
	upc_barrier;
	/* Time the construction of the hashtable */
	if (MYTHREAD == 0) {
	
		start = UPC_TICKS_NOW();

	}
#endif
	/* Build UFX hashtable */
   
#ifdef LSH
   int64_t oracleEntries;
   if (sscanf (string_size, "%lld", &oracleEntries)!=1) { printf ("ERROR: Size of oracle entries is not int64_t\n"); }
   loadOracleTable(oracleEntries, oracle_file , &oracle_table);
   upc_barrier;
   
#ifdef DEBUG_ORACLE
   if (MYTHREAD == 166) {
      
      FILE *printedOracle = fopen_chk("printedOracle", "w+");
      int64_t ite;
      char *localTable = (char*) malloc_chk(oracleEntries * sizeof(char));
      
      for (ite = 0; ite < oracleEntries; ite++) {
         
         int64_t original_hashval = ite;
         int64_t chunk1 = original_hashval / cores_per_node;
         int64_t offset1 = original_hashval % cores_per_node;
         int64_t nid = MYTHREAD / cores_per_node;
         int64_t nodes = (int64_t) ceil((double) THREADS / (double)cores_per_node);
         int64_t loc = nid * cores_per_node + chunk1 * THREADS + offset1;
         localTable[ite] = oracle_table[loc];
      
      }
      
      fwrite(localTable, sizeof(oracle_t), oracleEntries, printedOracle);
      fclose(printedOracle);
   
   }
   
   upc_barrier;
#endif

#endif
  
   
#ifndef BIN_INPUT
   /* FIXME: This works but is unnecessarily slow -- loook how to fix the siz ein the buildUFX funciton as did for the binary input case */
	size = load_factor * getNumKmersInUFX(input_UFX_name);
	if (MYTHREAD == 0) {
		int64_t minMemory = 12*(sizeof(list_t) + sizeof(int64_t)) * size / 10 / 1024/ 1024/ THREADS;
		printf("Minimum required shared memory: %lld MB. (%lld ufx kmers) If memory runs out re-run with more -shared-heap\n", minMemory, size);
	}
   fd = fopen_chk(input_UFX_name, "r");
   dist_hashtable = buildUFXhash(size, fd, &memory_heap, chunk_size);
   if (MYTHREAD == 0) printf("Number of buckets is %lld\n", dist_hashtable->size);
#endif
   

#ifdef BIN_INPUT
   int from_shm = 0;
   if (strstr(base_dir, "/dev/shm"))
       from_shm = 1;

   FILE *fdInput;
   int64_t myshare;
   int dsize;
   fdInput = UFXInitOpen(input_UFX_name, &dsize, &myshare, THREADS, MYTHREAD, &size, from_shm);
   if (MYTHREAD == 0) {
      int64_t minMemory = 12*(sizeof(list_t) + sizeof(int64_t)) * size / 10 / 1024/ 1024/ THREADS;
      printf("Minimum required shared memory: %ld MB. (%ld ufx kmers) If memory runs out re-run with more -shared-heap\n", minMemory, size);
   }
   dist_hashtable = buildUFXhash(size, fdInput, &memory_heap, myshare, dsize, dmin, chunk_size, load_factor);
#endif
   
   
   
	
#ifdef PROFILE
	/* Time the construction of the hashtable */
	if (MYTHREAD == 0) {

		end = UPC_TICKS_NOW();

		con_time = UPC_TICKS_TO_SECS(end-start);
		
		printf("\n\n*********** OVERALL TIME BREAKDOWN ***************\n\n");

		printf("\nTime for constructing UFX hash table is : %f seconds\n", con_time);
	}
#endif
	
	/* UU-mer graph traversal */

	/* initialze global shared variables */
	if (MYTHREAD == 0) {
		contig_id = 0;
		contig_list = NULL;
		upc_fence;
	}
   
	{
        	char outputfile_filepath[MAX_FILE_PATH];
	        sprintf(outputfile_filepath, "%s/%s_%d.fasta", base_dir, output_name, MYTHREAD);
		char *outputfile_name = get_rank_path(outputfile_filepath, MYTHREAD);
		my_out_file = fopen_chk(outputfile_name,"w");
	}
   
#ifdef PROFILE
	upc_barrier;
	
	/* Time the UU-graph traversal */
	if (MYTHREAD == 0) {
		start = UPC_TICKS_NOW();
	}
#endif

	UU_graph_traversal(dist_hashtable, my_out_file, minimum_contig_length);
	
#ifdef PROFILE
	upc_barrier;
	/* Time the UU-graph traversal */
	if (MYTHREAD == 0) {

		end = UPC_TICKS_NOW();

		trav_time = UPC_TICKS_TO_SECS(end-start);
		
		printf("\nTime for UU-graph traversal is : %f seconds\n", trav_time);
		printf("\nTotal time is :  %f seconds\n", trav_time + con_time);
		printf("\nTotal contigs stored: %ld\n", contig_id);
	}
#endif
   
#ifdef DETAILED_BUILD_PROFILE
   /* Measure load imbalance in blastmap */
   upc_barrier;
   ioTimes = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   storeTimes = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   calcTimes = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   setupTimes = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   upc_barrier;
   ioTimes[MYTHREAD] = fileIOTime;
   setupTimes[MYTHREAD] = setupTime;
   storeTimes[MYTHREAD] = storeTime;
   calcTimes[MYTHREAD] = cardCalcTime;
   upc_barrier;
   
   if (MYTHREAD == 0) {
      min_IO = ioTimes[0];
      max_IO = ioTimes[0];
      min_setup = setupTimes[0];
      max_setup = setupTimes[0];
      min_store = storeTimes[0];
      max_store = storeTimes[0];
      min_calc = calcTimes[0];
      max_calc = calcTimes[0];
      
      
      for (i=0; i<THREADS; i++) {
         sum_ios += ioTimes[i];
         sum_stores += storeTimes[i];
         sum_setups += setupTimes[i];
         sum_calcs += calcTimes[i];
         
         if ( ioTimes[i] < min_IO) {
            min_IO = ioTimes[i];
         }
         
         if ( ioTimes[i] > max_IO) {
            max_IO = ioTimes[i];
         }
         
         if ( storeTimes[i] < min_store) {
            min_store = storeTimes[i];
         }
         
         if ( storeTimes[i] > max_store) {
            max_store = storeTimes[i];
         }
         
         if ( calcTimes[i] < min_calc) {
            min_calc = calcTimes[i];
         }
         
         if ( calcTimes[i] > max_calc) {
            max_calc = calcTimes[i];
         }
         
         if ( setupTimes[i] < min_setup) {
            min_setup = setupTimes[i];
         }
         
         if ( setupTimes[i] > max_setup) {
            max_setup = setupTimes[i];
         }
         
      }
      printf("\n******** CONSTRUCTION DETAILED STATISTICS ********\n");
      printf("IOs: Avg %f\tMin %f\tMax %f\n",sum_ios/THREADS, min_IO, max_IO);
      printf("Calculating cardinalities: Avg %f\tMin %f\tMax %f\n",sum_calcs/THREADS, min_calc, max_calc);
      printf("Seting up: Avg %f\tMin %f\tMax %f\n",sum_setups/THREADS, min_setup, max_setup);
      printf("Storing k-mers: Avg %f\tMin %f\tMax %f\n",sum_stores/THREADS, min_store, max_store);
#ifdef LHS_PERF
      printf("Percentage of node-local lookups: %f %%\n", (100.0*success)/(1.0 * lookups));
#endif
      
   }
   
   upc_barrier;

#endif

	upc_barrier;
	fclose(my_out_file);
	
#ifdef PROFILE
#ifdef WORKLOAD_PROFILING
	if (MYTHREAD == 0) {
		printf("\n------------------------- STATS for UU graph traversal ----------------------------\n");
	}
	upc_barrier;
	fclose(my_out_file);
	printf("Thread %d added %lld bases in total\n", MYTHREAD, bases_added);
#endif
#endif

#ifdef UU_TRAV_PROFILE
	printf("Thread %d spent %f seconds in UU-graph traversal, %f seconds in walking (%f %%)\n", MYTHREAD, UU_time, walking_time, walking_time/UU_time*100.0);
#endif

   upc_barrier;

   if (mylog != NULL)
	   fclose(mylog);
   
   /* Print some metafiles used in downstream steps */
   if ( MYTHREAD == 0 ) {
      char countFileName[MAX_FILE_PATH];
      sprintf(countFileName,"n%s.txt", output_name);
      FILE *countFD = fopen_rank_path(countFileName, "w+", -1);
      fprintf(countFD, "%ld\n", contig_id);
      fclose(countFD);
   }
   
   {
      char mycountFileName[MAX_FILE_PATH];
      sprintf(mycountFileName,"%s/my%s_%d.txt", base_dir, output_name, MYTHREAD);
      FILE *mycountFD = fopen_rank_path(mycountFileName, "w+", MYTHREAD);
      fprintf(mycountFD, "%d\n", myContigs);
      fclose(mycountFD);
   }
   
   upc_barrier;

   if (!MYTHREAD)
	   printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
			  ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));

	return 0;
}
