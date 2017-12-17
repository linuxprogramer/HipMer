#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "../common/optlist.h"
#include <upc.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include "fasta.c"
#include <errno.h>
#include <upc_tick.h>
#include <libgen.h>

#include "../common/common.h"
#include "../common/upc_compatibility.h"
#include "../contigs/kmer_handling.h"
#include "kmer_hash2.h"
#include "contigDB.h"
#include "buildSeedIndex.h"

#ifdef PROFILE
double setup_cache = 0.0;
double readmap_time = 0.0;
double fetch_contigs_time = 0.0;
double hash_table_lookups = 0.0;
double input_IO = 0.0;
double output_IO = 0.0;
#endif

int64_t local_allocs;
int64_t cache_hits;
int64_t cache_queries;
int64_t seedsExtracted = 0;
int cores_per_node = 24;

#define NOT_EXCLUDE
#ifdef NOT_EXCLUDE
#include "aligningPhase.h"
#endif

shared int64_t globalSeedCardinality = 0;
shared int64_t globalSubcontigCardinality = 0;
shared int64_t globalAlignmentsFound = 0;
shared int64_t globalReadsProcessed = 0;
shared int64_t globalReadsMapped= 0;
shared int64_t globalLooks= 0;
shared int64_t globalSeedsExtracted=0;

int main(int argc, char **argv) {
   upc_tick_t start_time = upc_ticks_now();
   int readNumber;
   int i, Fsubcontigs;
   int64_t size;
   hash_table_t *dist_hashtable;
   memory_heap_t memory_heap;
   FILE *fd, *my_out_file, *blastresFD1, *blastresFD2, *cacheFD;
   FILE *reads_fd, *fdReads;
   char readFileName[MAX_FILE_PATH];
   char outputBlastmap1[MAX_FILE_PATH], outputBlastmap2[MAX_FILE_PATH];
   shared[1] double *fetch_contigs_times;
   shared[1] double *hash_table_lookups_times;
   shared[1] double *readmap_times;
   shared[1] double *outputIO_times;
   double min_IO, min_fetch, min_lookup, min_readmap;
   double max_IO, max_lookup, max_fetch, max_readmap;
   double sum_lookups=0.0;
   double sum_outputIO=0.0;
   double sum_fetchcontigs=0.0;
   double sum_readmaps=0.0;
   int64_t nSlots;
   int64_t mySeedCardinality, partial_result, FseedCardinality;
   int64_t localSubContigCardinality, IDoffset;
   
   /* Use getopt() to read arguments */
   //extern char *optarg;
   //extern int optind;
   shared[1] list_t *cacheTable;
   shared[1] int64_t *cacheFlags;

   int c, filesPerPair;
   const char *base_dir = ".";
   char *libName, *read_files_name, *cache_size_in_MB, *minimum_contig_length, *UFX_boundaries, *dmin, *nContigs, *chunksize, *kmer_cache_in_MB, *exp_factor_string, *contigFileName;
   double exp_factor;
   int read_length = 0;
   int scEffLength = 0;
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "l:P:r:c:C:m:b:d:n:x:k:e:j:O:B:N:L:");
   filesPerPair = 1; // default interleaved
   
   int outputSize;
   libName = read_files_name = NULL;
   
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'l':
            libName = thisOpt->argument; 
            break;
         case 'P':
            filesPerPair = atoi(thisOpt->argument);
            break;
         case 'e':
            scEffLength = atoi(thisOpt->argument);
            break;
         case 'r':
            read_files_name = thisOpt->argument;
            break;
         case 'C':
            cache_size_in_MB = thisOpt->argument;
            break;
         case 'c':
            contigFileName = thisOpt->argument;
            break;
         case 'm':
            minimum_contig_length = thisOpt->argument;
            break;
         case 'b':
            UFX_boundaries = thisOpt->argument;
            break;
         case 'd':
            dmin = thisOpt->argument;
            break;
         case 'n':
            nContigs = thisOpt->argument;
            break;
         case 'x':
            chunksize = thisOpt->argument;
            break;
         case 'k':
            kmer_cache_in_MB = thisOpt->argument;
            break;
         case 'j':
            exp_factor_string = thisOpt->argument;
            break;
         case 'O':
            outputSize = atoi(thisOpt->argument);
            break;
         case 'B':
             base_dir = thisOpt->argument;
             break;
         case 'N':
			 cores_per_node = atoi(thisOpt->argument);
			 break;
         case 'L':
             read_length = atoi(thisOpt->argument);
             break;
         default:
            break;
      }
      
      free(thisOpt);
   }

   if (!MYTHREAD)
	   printf("Using %d cores per node\n", cores_per_node);

   if (!read_length) 
       scEffLength = 0;
   
    UPC_TICK_T start, end, start_timer, end_timer;
    double con_time, trav_time, blastmap_time;
    char referenceFileName[MAX_FILE_PATH];

   if (libName == NULL || read_files_name  == NULL) 
       DIE("you must specifiy -l and -r (%s %s)!\n", 
           !libName ? "NULL" : libName, !read_files_name ? "NULL" : read_files_name);
   /* TODO check other necessary parameters */
   
   sprintf(referenceFileName,"%s/%s_%d.fasta", base_dir, contigFileName, MYTHREAD);
   get_rank_path(referenceFileName, MYTHREAD);
   
   if (strncmp(base_dir, "/dev/shm", 8) == 0) {
        // FIXME: crude hack for removing ufx files from shm. This probably shouldn't be done here, but 
        // somewhere else in the pipeline
        char ufx_name[MAX_FILE_PATH];
        sprintf(ufx_name, "/dev/shm/ALL_INPUTS.fofn.ufx.bin%d", MYTHREAD);
        unlink(ufx_name);
        //printf("thread %d unlinking %s\n", MYTHREAD, ufx_name);
        sprintf(ufx_name, "/dev/shm/ALL_INPUTS.fofn.ufx.bin%d.entries", MYTHREAD);
        unlink(ufx_name);
        //printf("thread %d unlinking %s\n", MYTHREAD, ufx_name);
   }
    
   char cache_statistics_name[MAX_FILE_PATH];
   sprintf(cache_statistics_name,"cache_statistics_%d", MYTHREAD);
   get_rank_path(cache_statistics_name, MYTHREAD);
   
   /************************************************/
   /* Estimate an upper bound for seed cardinality */
   /************************************************/
   exp_factor = atof(exp_factor_string);
   
#ifdef PROFILE
	upc_barrier;
	if (MYTHREAD == 0) {
      start = UPC_TICKS_NOW();
	}
#endif

   mySeedCardinality = seedAndContigCardinality(referenceFileName, &localSubContigCardinality, scEffLength);
   partial_result = bupc_atomicI64_fetchadd_strict(&globalSeedCardinality, mySeedCardinality);
   IDoffset = bupc_atomicI64_fetchadd_strict(&globalSubcontigCardinality, localSubContigCardinality);

   upc_barrier;
   FseedCardinality = globalSeedCardinality;
   Fsubcontigs = globalSubcontigCardinality;
   upc_fence;
   
#ifdef PROFILE
   upc_barrier;
	if (MYTHREAD == 0) {
      end = UPC_TICKS_NOW();
		con_time = (UPC_TICKS_TO_SECS(end-start));
		printf("\nTime for calculating seed cardinality is : %f seconds\n", con_time);
	}
#endif
   
   /************************/
	/*   Build seed index   */
   /************************/

#ifdef PROFILE
	upc_barrier;
	if (MYTHREAD == 0) {
      start = UPC_TICKS_NOW();
	}
#endif
   
    dist_hashtable = buildSeedIndex(FseedCardinality, exp_factor, referenceFileName, &memory_heap, scEffLength, IDoffset, localSubContigCardinality, read_length);
   
#ifdef PROFILE
   upc_barrier;
	if (MYTHREAD == 0) {
      end = UPC_TICKS_NOW();
		con_time = (UPC_TICKS_TO_SECS(end-start));
		printf("\nTime for building the seed index is : %f seconds\n", con_time);
	}
#endif
      
   /*****************************/
	/* Parallel sequence aligner */
   /*****************************/
#ifdef NOT_EXCLUDE
   
#ifdef PROFILE
   upc_barrier;
   if (MYTHREAD == 0) {
      printf("\n\n***************** Starting aligning *******************\n\n");
		start = UPC_TICKS_NOW();
   }
#endif

   
#ifdef PROFILE
   start_timer = UPC_TICKS_NOW();
#endif
   /* Build cache for k-mer seeds */
#ifdef USE_SWCACHE
   upc_barrier;
   int64_t expanded_size = ((int64_t) atoi(kmer_cache_in_MB)) * ((int64_t) (1024 * 1024));
   nSlots = createKmerCache(&cacheTable, expanded_size, &cacheFlags );
   upc_barrier;
#endif

   shared[1] contigDataPtr *cacheTableContig = NULL;
   shared int64_t *cachePtr = NULL;
#ifdef USE_SWCACHE
   int64_t expanded_size2 =((int64_t) atoi(cache_size_in_MB)) * ((int64_t) 1024 * 1024);
   createCache(Fsubcontigs, &cachePtr, &cacheTableContig, expanded_size2);
#endif
   
#ifdef PROFILE
   end_timer = UPC_TICKS_NOW();
   setup_cache += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
   
   int64_t readsMapped = 0, readsProcessed = 0, alignmentsFound = 0;
   
   /* delete any preexisting merAlignerOutput file */

   sprintf(outputBlastmap1, "%s/%s-merAlignerOutput_%d_Read1", base_dir, libName, MYTHREAD );
   get_rank_path(outputBlastmap1, MYTHREAD);
   blastresFD1 = fopen_chk(outputBlastmap1,"w");
   if (filesPerPair > 0) {
       sprintf(outputBlastmap2, "%s/%s-merAlignerOutput_%d_Read2", base_dir, libName, MYTHREAD );
       get_rank_path(outputBlastmap2, MYTHREAD);
       blastresFD2 = fopen_chk(outputBlastmap2,"w");
   } else {
       blastresFD2 = NULL;
   }

   readNumber=0;
   fdReads = fopen_chk(read_files_name, "r");
   while ( fgets(readFileName, MAX_FILE_PATH, fdReads) != NULL ) {
      readFileName[strlen(readFileName)-1] = '\0';
      if ( readFileName[0] == '\0' ) break ;
      if (filesPerPair == 2) {
        /* switch read number */
        readNumber = (readNumber+1) % 2; 
      }

      alignmentsFound += parallelAligner(dist_hashtable , readFileName, blastresFD1, blastresFD2, filesPerPair, readNumber, 0, 0, Fsubcontigs, atoi(cache_size_in_MB), atoi(minimum_contig_length), atoi(chunksize), nSlots, cacheTable, cacheFlags, cachePtr, cacheTableContig, &readsMapped, &readsProcessed, base_dir, libName);
   }

   fclose(blastresFD1);
   if (filesPerPair > 0) {
     fclose(blastresFD2);
   } 
   fclose(fdReads);

#ifdef PROFILE
   if (MYTHREAD == 0) {
      end = UPC_TICKS_NOW();
      blastmap_time = (UPC_TICKS_TO_SECS(end-start));
      printf("Software contig cache ENABLED with size %f GB/node\n", atoi(cache_size_in_MB)/1024.0);
      printf("Software k-mer cache ENABLED with size %f GB/node\n", atoi(kmer_cache_in_MB)/1024.0);
      printf("\nParallel merAligner took :  %f seconds\n", blastmap_time);
   }
#endif

   upc_barrier;
   
   /* Measure load imbalance in aligning phase */
   readmap_times = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   outputIO_times = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   hash_table_lookups_times = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   fetch_contigs_times = (shared[1] double*) upc_all_alloc(THREADS, sizeof(double));
   if (readmap_times == NULL || outputIO_times == NULL || hash_table_lookups_times == NULL || fetch_contigs_times == NULL) 
       DIE("Could not allocate timing datastructures\n");
   upc_barrier;
   readmap_times[MYTHREAD] = readmap_time;
   outputIO_times[MYTHREAD] = output_IO;
   hash_table_lookups_times[MYTHREAD] = hash_table_lookups;
   fetch_contigs_times[MYTHREAD] = fetch_contigs_time;
   
   partial_result = bupc_atomicI64_fetchadd_strict(&globalAlignmentsFound, alignmentsFound);
   partial_result = bupc_atomicI64_fetchadd_strict(&globalReadsMapped, readsMapped);
   partial_result = bupc_atomicI64_fetchadd_strict(&globalReadsProcessed, readsProcessed);
   partial_result = bupc_atomicI64_fetchadd_strict(&globalLooks, cache_queries);
   partial_result = bupc_atomicI64_fetchadd_strict(&globalSeedsExtracted, seedsExtracted);


   upc_barrier;
   
   if (MYTHREAD == 0) {
      min_fetch = fetch_contigs_times[0];
      max_fetch = fetch_contigs_times[0];
      min_readmap = readmap_times[0];
      max_readmap = readmap_times[0];
      min_IO = outputIO_times[0];
      max_IO = outputIO_times[0];
      min_lookup = hash_table_lookups_times[0];
      max_lookup = hash_table_lookups_times[0];
      
      
      for (i=0; i<THREADS; i++) {
         sum_readmaps += readmap_times[i];
         sum_outputIO += outputIO_times[i];
         sum_lookups += hash_table_lookups_times[i];
         sum_fetchcontigs += fetch_contigs_times[i];
         
         if ( readmap_times[i] < min_readmap) {
            min_readmap = readmap_times[i];
         }
         
         if ( readmap_times[i] > max_readmap) {
            max_readmap = readmap_times[i];
         }
         
         if ( outputIO_times[i] < min_IO) {
            min_IO = outputIO_times[i];
         }
         
         if ( outputIO_times[i] > max_IO) {
            max_IO = outputIO_times[i];
         }
         
         if ( hash_table_lookups_times[i] < min_lookup) {
            min_lookup = hash_table_lookups_times[i];
         }
         
         if ( hash_table_lookups_times[i] > max_lookup) {
            max_lookup = hash_table_lookups_times[i];
         }
         
         if ( fetch_contigs_times[i] < min_fetch) {
            min_fetch = fetch_contigs_times[i];
         }
         
         if ( fetch_contigs_times[i] > max_fetch) {
            max_fetch = fetch_contigs_times[i];
         }
         
      }
      printf("\n*******  AGGREGATE STATISTICS  *******\n");
      printf("Readmaps: Avg %f\tMin %f\tMax %f\n",sum_readmaps/THREADS, min_readmap, max_readmap);
      printf("Look-ups: Avg %f\tMin %f\tMax %f\n",sum_lookups/THREADS, min_lookup, max_lookup);
      printf("Fetching contigs: Avg %f\tMin %f\tMax %f\n",sum_fetchcontigs/THREADS, min_fetch, max_fetch);
      printf("Output IO: Avg %f\tMin %f\tMax %f\n",sum_outputIO/THREADS, min_IO, max_IO);
      
      printf("\n**************  ALIGNMENTS STATISTICS  ********************\n");
      printf("Total reads processed: %ld\n",globalReadsProcessed);
      printf("Total reads mapped: %ld (%f %%)\n",globalReadsMapped, ((double)globalReadsMapped/(double)globalReadsProcessed)*100.0);
      printf("Total alignments found: %ld\n",globalAlignmentsFound);
      printf("Global seed lookups: %ld\n", globalLooks);
      printf("Global seed cardinality is: %ld\n", globalSeedCardinality);
      printf("Global subcontig cardinality is: %ld\n", globalSubcontigCardinality);
      int64_t tmpres = 0;
      for (i=0; i<THREADS; i++) {
         tmpres += memory_heap.heap_indices[i];
      }
      
      {
          char countFileName[MAX_FILE_PATH];
          sprintf(countFileName,"%s-nTotalAlignments.txt", libName);
          char *fname = get_rank_path(countFileName, -1);
          FILE *countFD = fopen_chk(fname, "w" );
          fprintf(countFD, "%ld\n", globalAlignmentsFound);
          fclose(countFD);
      }
      
      printf("Entries in distributed seed index are: %ld\n", tmpres);
      printf("Number of actual seeds extracted is: %ld\n", globalSeedsExtracted);
      
   }
   
   upc_barrier;
   //cacheFD = fopen_chk(cache_statistics_name,"w+");
   //fprintf(cacheFD, "THREAD %d did %lld k-mer queries and cache hits were %lld  -- hit ratio is %f\n", MYTHREAD, cache_queries, cache_hits, (double)cache_hits/(double) cache_queries);
   //fclose(cacheFD);
   
   /* Dump files with detailed statistics */
   
#endif // ifdef NOT_EXCLUDE


   char detailed_statistics_name[255];
   FILE *detailedFD;
	sprintf(detailed_statistics_name,"detailed_statistics_%d", MYTHREAD);
   //detailedFD = fopen_chk(detailed_statistics_name,"w+");
   //fprintf(detailedFD, "THREAD %d spent:\nReadmap time = %f seconds\n=========================\nSeed lookups = %f seconds\nFetching contigs = %f seconds\nComputation time = %f seconds\nOutput I/O = %f seconds\n", MYTHREAD, readmap_time, hash_table_lookups, fetch_contigs_time, readmap_time - hash_table_lookups - fetch_contigs_time - output_IO, output_IO );
   //fprintf(detailedFD, "Input I/O = %f seconds\n", input_IO);
   //fclose(detailedFD);
   upc_barrier;
   
	if (!MYTHREAD)
		printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
			   ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
   return 0;
}
