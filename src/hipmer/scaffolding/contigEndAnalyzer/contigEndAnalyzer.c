#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <assert.h>
#include <libgen.h>
#include "meraculous.h"
#include "../../common/optlist.h"
#include "../../common/common.h"
#include "../../common/kseq.h"
#include "../../common/Buffer.h"
#include "../../contigs/kmer_hash.h"
#include "../../kmercount/readufx.h"

// these are used in buildUFXHashBinary.h and so have to be defined before it is included
double fileIOTime = 0.0;
double cardCalcTime = 0.0;
double setupTime = 0.0;
double storeTime = 0.0;

#include "../../contigs/buildUFXhashBinary.h"

#define MAX_CONTIG_SIZE (1<<20)

int main(int argc, char **argv) {
    upc_tick_t start_time = upc_ticks_now();
   
   /* Use getopt() to read arguments */
   extern char *optarg;
	extern int optind;
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "i:c:k:s:d:B:");
   
   char *input_UFX_name, *contig_file_name;
   int kmerLength;
   int nContigs, dmin, chunk_size=1;

   char *base_dir = ".";

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
         case 'd':
            dmin = atoi(thisOpt->argument);
            break;
         case 'c':
            contig_file_name = thisOpt->argument;
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
   sprintf(logfile_name,"contigEndAnalyser_%s_%d.log", contig_file_name, MYTHREAD);
   
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
   
   /* Build hash table using UFX file that also contains the FX kmers */
   
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

   /* Read contigs and find mean depth for each one */

   gzFile contigFile;
   FILE *myOutputFile;
   int contigID, is_least;
   char *firstKmer, *lastKmer, *prevKmer, *nextKmer;
   list_t copy;
   shared[] list_t *lookup_res = NULL;
   int nMers = 0, contigLen, i;
   char prevBase, nextBase, prevCodeL, prevCodeR, nextCodeL, nextCodeR;
   
   firstKmer = (char *) malloc_chk((kmerLength+1) * sizeof(char));
   firstKmer[kmerLength] = '\0';
   lastKmer = (char *) malloc_chk((kmerLength+1) * sizeof(char));
   lastKmer[kmerLength] = '\0';
   prevKmer = (char *) malloc_chk((kmerLength+1) * sizeof(char));
   prevKmer[kmerLength] = '\0';
   nextKmer = (char *) malloc_chk((kmerLength+1) * sizeof(char));
   nextKmer[kmerLength] = '\0';
   
   {
       char my_contig_file_name[MAX_FILE_PATH];
       sprintf(my_contig_file_name,"%s/%s_%d.fasta", base_dir, contig_file_name, MYTHREAD);
       get_rank_path(my_contig_file_name, MYTHREAD);
       contigFile = gzopen(my_contig_file_name, "r");
       if (!contigFile) { DIE("Could not open %s\n", my_contig_file_name); }
   }
   {
       char my_output_file_name[MAX_FILE_PATH];
       sprintf(my_output_file_name,"%s/%s_%d.cea", base_dir, contig_file_name, MYTHREAD);
       get_rank_path(my_output_file_name, MYTHREAD);
       myOutputFile = fopen_chk(my_output_file_name, "w");
   }
   
   kseq_t *ks = kseq_init(contigFile);
   while ( kseq_read(ks) >= 0 ) {
      
      /* Read a contig and its length */
      contigID = atoi(ks->name.s+7);

      contigLen = ks->seq.l;
      assert(contigLen >= kmerLength);
      
      /* Extract the first and prev kmer of the contig */
      memcpy(firstKmer, ks->seq.s, kmerLength * sizeof(char));
      lookup_res = lookup_least_kmer_and_copy(dist_hashtable, firstKmer, &copy, &is_least);
      if (is_least) {
         prevBase = copy.left_ext;
      } else {
         prevBase = reverseComplementBaseExt(copy.right_ext);
      }
      prevKmer[0] = prevBase;

      memcpy(&prevKmer[1], firstKmer, (kmerLength-1) * sizeof(char));
      lookup_res = lookup_least_kmer_and_copy(dist_hashtable, prevKmer, &copy, &is_least);
      if (lookup_res != NULL) {
         if (is_least) {
            prevCodeL = copy.left_ext;
            prevCodeR = copy.right_ext;
         } else {
            prevCodeL = reverseComplementBaseExt(copy.right_ext);
            prevCodeR = reverseComplementBaseExt(copy.left_ext);
         }
      } else {
         prevCodeL = '0';
         prevCodeR = '0';
      }
      
      /* Extract the last and next kmer of the contig */
      memcpy(lastKmer, ks->seq.s + contigLen-kmerLength, kmerLength * sizeof(char));
      lookup_res = lookup_least_kmer_and_copy(dist_hashtable, lastKmer, &copy, &is_least);
      if (is_least) {
         nextBase = copy.right_ext;
      } else {
         nextBase = reverseComplementBaseExt(copy.left_ext);
      }
      nextKmer[kmerLength-1] = nextBase;
      memcpy(nextKmer, &lastKmer[1], (kmerLength-1) * sizeof(char));
      lookup_res = lookup_least_kmer_and_copy(dist_hashtable, nextKmer, &copy, &is_least);
      if (lookup_res != NULL) {
         if (is_least) {
            nextCodeL = copy.left_ext;
            nextCodeR = copy.right_ext;
         } else {
            nextCodeR = reverseComplementBaseExt(copy.left_ext);
            nextCodeL = reverseComplementBaseExt(copy.right_ext);
         }
      }  else {
         nextCodeL = '0';
         nextCodeR = '0';
      }

      /* Print to the output file */
      fprintf(myOutputFile,"Contig%d\t[%c%c]\t(%c)\t%s\t%d\t%s\t(%c)\t[%c%c]\n", contigID, prevCodeL, prevCodeR, prevBase, firstKmer, contigLen, lastKmer, nextBase, nextCodeL, nextCodeR );
   }
   
   kseq_destroy(ks);
   gzclose(contigFile);
   fclose(myOutputFile);
   
   upc_barrier;
   if (MYTHREAD == 0) {
      end = UPC_TICKS_NOW();
      depth_time = UPC_TICKS_TO_SECS(end-start);
      printf("\nTime for analyzing the contig ends : %f seconds\n", depth_time);
   }
   
   if (!MYTHREAD)
       printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
             ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));

   return 0;
}
