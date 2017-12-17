#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <assert.h>
#include <libgen.h>
#include "../../common/optlist.h"
#include "../../common/common.h"
#include "../../common/Buffer.h"
#include "../../contigs/kmer_hash.h"
#include "../../kmercount/readufx.h"
#include "../../common/kseq.h"
#include "meraculous.h"

double fileIOTime = 0.0;
double cardCalcTime = 0.0;
double setupTime = 0.0;
double storeTime = 0.0;

#include "../../contigs/buildUFXhashBinary.h"

#define MAX_CONTIG_SIZE 900000

int main(int argc, char **argv) {
    upc_tick_t start_time = upc_ticks_now();
   
   /* Use getopt() to read arguments */
   extern char *optarg;
   extern int optind;
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "i:c:k:d:s:B:");
   
   char *input_UFX_name, *contig_file_name;
   int kmerLength, dmin, chunk_size;

   const char *base_dir = ".";

   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'i':
            input_UFX_name = thisOpt->argument;
            break;
         case 'k':
            kmerLength = atoi(thisOpt->argument);
            break;
         case 'c':
            contig_file_name = thisOpt->argument;
            break;
         case 'd':
            dmin = atoi(thisOpt->argument);
            break;
         case 's':
            chunk_size = atoi(thisOpt->argument);
            break;
         case 'B':
             base_dir = thisOpt->argument;
             break;
         default:
            break;
      }
      
      free(thisOpt);
   }
   
   hash_table_t *dist_hashtable;
   UPC_TICK_T start, end;
   
   if (MYTHREAD == 0) {
#ifndef NO_PAD
       printf("Struct size is %lu, with %lu padding, %lu shared[] ptr\n", 
              (unsigned long)sizeof(list_t), (unsigned long)UPC_PADDING, (unsigned long)sizeof(shared void*));
#else
      printf("Struct size is %lu, no padding, %lu shared[] ptr\n", 
             (unsigned long)sizeof(list_t), (unsigned long)sizeof(shared void*));
#endif
   }
   double con_time, depth_time;
   
#if VERBOSE > 0
   char logfile_name[MAX_FILE_PATH];
   sprintf(logfile_name,"contigMerDepth_%s_%d.log", contig_file_name, MYTHREAD);
   
   mylog = fopen_rank_path(logfile_name, "w", MYTHREAD);
#else
   mylog = NULL;
#endif
   
   
#ifdef PROFILE
   upc_barrier;
   /* Time the construction of the hashtable */
   if (MYTHREAD == 0) {
      start = UPC_TICKS_NOW();
   }
#endif

   /* Build hash table using UFX file that also contains k-mer depths */
   
#ifdef BIN_INPUT
   FILE *fdInput;
   int64_t myshare;
   int dsize;
   int64_t size;
   memory_heap_t memory_heap;
   int from_shm = 0;
   if (strstr(base_dir, "/dev/shm"))
       from_shm = 1;
   fdInput = UFXInitOpen(input_UFX_name, &dsize, &myshare, THREADS, MYTHREAD, &size, from_shm);
   if (MYTHREAD == 0) {
      int64_t minMemory = 12*(sizeof(list_t) + sizeof(int64_t)) * size / 10 / 1024/ 1024/ THREADS;
      printf("Minimum required shared memory: %ld MB. (%ld ufx kmers) If memory runs out re-run with more -shared-heap\n", minMemory, size);
   }
   dist_hashtable = buildUFXhash(size, fdInput, &memory_heap, myshare, dsize, dmin, chunk_size, 1);
#endif
   
   
#ifdef PROFILE
   upc_barrier;
   /* Time the construction of the hashtable */
   if (MYTHREAD == 0) {
      end = UPC_TICKS_NOW();
      con_time = UPC_TICKS_TO_SECS(end-start);
      printf("\n\n*********** OVERALL TIME BREAKDOWN ***************\n\n");
      printf("\nTime for constructing UFX hash table is : %f seconds\n", con_time);
      start = UPC_TICKS_NOW();
   }
#endif
   
   /* Calculating the depth of the contigs */
   
   /* Read contigs and find mean depth for each one */
   char my_contig_file_name[MAX_FILE_PATH];
   char my_output_file_name[MAX_FILE_PATH];

   FILE *myOutputFile;
   gzFile contigFile;
   int contigID, is_least;
   char *curKmer;
   list_t copy;
   shared[] list_t *lookup_res = NULL;
   int64_t runningDepth = 0;
   int nMers = 0, contigLen, i;
   double mean = 0.0;
   
   curKmer = (char *) malloc_chk((kmerLength+1) * sizeof(char));
   curKmer[kmerLength] = '\0';
   
   sprintf(my_contig_file_name,"%s/%s_%d.fasta", base_dir, contig_file_name, MYTHREAD);
   get_rank_path(my_contig_file_name, MYTHREAD);
   contigFile = gzopen(my_contig_file_name, "r");
   if (!contigFile) { DIE("Could not open %s\n", my_contig_file_name); }
   kseq_t *ks = kseq_init(contigFile);

   sprintf(my_output_file_name,"%s/merDepth_%s_%d.txt", base_dir, contig_file_name, MYTHREAD);
   myOutputFile = fopen_rank_path(my_output_file_name, "w", MYTHREAD);

   while ( kseq_read(ks) >= 0 ) {
      /* Read a contig and its length */
      contigID = atoi(ks->name.s+7);
      contigLen = ks->seq.l;
      runningDepth = 0;
      /* For each kmer in the contig extract the depth and add to the running sum */
      for (i=0; i<=contigLen-kmerLength; i++) {
         memcpy(curKmer, ks->seq.s + i, kmerLength * sizeof(char));
         lookup_res = lookup_least_kmer_and_copy(dist_hashtable, curKmer, &copy, &is_least);
         runningDepth += copy.count;
      }
      
      nMers = contigLen-kmerLength+1;
      mean = (1.0 * runningDepth)/(1.0 * nMers);
      fprintf(myOutputFile, "Contig%d\t%d\t%f\n", contigID, nMers, mean);
   }
   
   kseq_destroy(ks);
   gzclose(contigFile);
   fclose(myOutputFile);
   
   upc_barrier;
   if (MYTHREAD == 0) {
      end = UPC_TICKS_NOW();
      depth_time = UPC_TICKS_TO_SECS(end-start);
      printf("\nTime for calculating the contig depths : %f seconds\n", depth_time);
   }
   
   if (!MYTHREAD)
       printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
              ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));

   return 0;
}
