#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <upc.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <libgen.h>

#include "../../common/optlist.h"
#include "../../common/common.h"
#include "../../common/upc_compatibility.h"
#include "bma_meta.h"
#include "gapSizeEstimatesUtils.h"
#include "linkHash.h"
int scaffoldSpans = 0;
#include "bmaToLinksUtils.h"

shared int64_t nTRUNC = 0;
shared int64_t nSMALL = 0;
shared int64_t nINORM = 0;
shared int64_t nINNIE = 0;
shared int64_t nSELFL = 0;
shared int64_t nEDIST = 0;
shared int64_t nREDUN = 0;
shared int64_t nACCPT = 0;
shared int64_t nISHIRT = 0;

shared int64_t tot_links = 0;
shared int64_t tot_nLinks = 0;
shared int64_t PRINTnLinks = 0;
shared int64_t tot_spans = 0;

#define LINK_LEN 200
#define COL0_LEN 20
#define COL1_LEN 100
#define COL2_LEN 100
#define COL3_LEN 300
#define COL4_LEN 100
#define COL5_LEN 300

int main(int argc, char **argv) {
   upc_tick_t start_time = upc_ticks_now();
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "b:m:D:P:n:L:l:B:");
   UPC_TICK_T start, end;
   char *cur_bmaFile;
   int linkUid = -1;
   char *bmaDataFileNames = NULL;
   FILE *cur_bma_fd;
   int merSize, minNetEndDistance = 0, previousMaxInsertSize = 0;
   int64_t i, j, total_entries, nLibs = 0;
   char *lineBuffers = (char*) calloc_chk(MY_LINESIZE * 2 + LINK_LEN + COL0_LEN + COL1_LEN + COL2_LEN + COL3_LEN + COL4_LEN + COL5_LEN, 1);
   char *line = lineBuffers;
   char *aux_line = line + MY_LINESIZE;
   char *link = aux_line + MY_LINESIZE;
   char *col0 = link + LINK_LEN;
   char *col1 = col0 + COL0_LEN;
   char *col2 = col1 + COL1_LEN;
   char *col3 = col2 + COL2_LEN;
   char *col4 = col3 + COL3_LEN;
   char *col5 = col4 + COL4_LEN;
   char personalBmaDataFile[MAX_FILE_PATH];
   char personalBmaDataFileOut[MAX_FILE_PATH];
   FILE *outFD;
   bma_info bmas[MAX_LIB_NUM];
   int64_t *local_index;
   link_heap_t link_heap;
   int insertSize = -1, insertStdDev = 0, my_count = 0;
   /*  Each line for a library: column0 --> spans_count | column1 --> endSeparation sum */
   int64_t spans[MAX_LIB_NUM][3];
   int64_t *splints;
   int cpy_previousMaxInsertSize = 0;
   int64_t datum_heap_pos, link_heap_pos, entries_in_my_stack;
   datum_t *datum_heap, *cur_data;
   linkList_t *linkList_heap, *cur_bucket;
   link_hash_table_t *local_hashtable;
   link_t *local_buffs;
   link_t *my_stack;
   int minimumGap, splintMaxDev, spanMaxZ, cur_e1, cur_e2, cur_type, cur_lib, cur_endSep, splintDev,my_checksum;
   int64_t  splintMaxFreq ,nSplints, nSpans, anomalousSpans, anomalousSplints;
   int spanAnomaly, splintGapEstimate, nLibSpans, readLength, netSeparation, l1, l2, spanGapUncertainty;
   double spanZ;
   double meanGapEstimate, meanOffset, sumOfWeights, spanGapUncertaintyDouble, libWeight;
   char cur_compressedTails;
   int tail1, tail2, my_anomalies;
   double spansGE[MAX_LIB_NUM], spanGapEstimate, gapEstimate;
   int64_t total;
   int64_t myPRINTlinks = 0;
   
   char my_nContigsFile[255];
   FILE *my_nC_fd;

   const char *base_dir = ".";

   UPC_TICK_T startT, endT;
   UPC_TICK_T start_io, end_io;
   double io_read_time = 0.0;
   
   char *libUid;
   
   splints = (int64_t*) malloc_chk(MAX_END_SEPARATION * sizeof(int64_t));
   /* Process the input arguments */
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'b':
            bmaDataFileNames = thisOpt->argument;
            break;
         case 'm':
            merSize = atoi(thisOpt->argument);
            break;
         case 'D':
            minNetEndDistance = atoi(thisOpt->argument);
            break;
         case 'P':
            previousMaxInsertSize = atoi(thisOpt->argument);
            cpy_previousMaxInsertSize = previousMaxInsertSize;
            break;
         case 'n':
            /* Total number of entries in bma files */
            total_entries = __atoi64(thisOpt->argument);
            break;
         case 'L':
            linkUid = atoi(thisOpt->argument);
            break;
         case 'B':
             base_dir = thisOpt->argument;
             break;
         default:
            break;
      }
      free(thisOpt);
   }
   
   if (linkUid < 0 || total_entries <= 0 || !bmaDataFileNames) 
       DIE("Please enter -b bmaMetaFileNames -n (estimated)TotalEntries and" 
           "a L linkUid(round number)\nwhere bmaMetaFileNames is a comma-separated list\n");
   upc_barrier;
   start = UPC_TICKS_NOW();
   /* bmaToLinks code */

   char *next_name = NULL;
   char *cur_name = bmaDataFileNames;
   /* Count the number of libs and store info for each library in the mapInfo data structure */
   while (1) {
       assert(nLibs < MAX_LIB_NUM);
       char *next_name = strchr(cur_name, ',');
       if (next_name) 
           next_name[0] = 0;

       // one lib per file
       char name_path[MAX_FILE_PATH];
       strcpy(name_path, cur_name);
       FILE *bmaDataFile = fopen_rank_path(name_path, "r", -1);
       if (!fgets(line, MY_LINESIZE, bmaDataFile)) 
           DIE("Could not get the line in %s -> %s\n", cur_name, name_path);
       strcpy(aux_line, line);
       split_bmaLine(aux_line, &bmas[nLibs], nLibs);
       for (i=0 ; i < REASONS; i++)
           bmas[nLibs].libPairSummary[i] = 0;
#ifdef DEBUG
       fprintf(stderr, "Thread %d read file %s (%"PRId64 "): %s %d %d\n", 
                     MYTHREAD, name_path, nLibs, bmas[nLibs].bmaFile, bmas[nLibs].insertSize, bmas[nLibs].stdDev);
#endif
       fclose(bmaDataFile);
       nLibs++;
       if (!next_name)
           break;
       cur_name = next_name + 1;
   }
   
   upc_barrier;
   
   startT = UPC_TICKS_NOW();
   
   serial_printf("Threads done reading BMA Metafiles, found %ld libs\n", nLibs);
   
   /* Over-allocating in heaps to avoid OOM errors in case of large load imbalance */
   total_entries = total_entries * EXP_FACTOR;
   
   /* Allocate local buffers and local shared heaps required for the links hashtable */
   allocate_local_buffs(&local_buffs, &local_index);
   create_link_heaps(total_entries, &link_heap);
   
   char *fgets_result;
   
   upc_barrier;
   serial_printf("Threads done allocating link heaps: %ld\n", total_entries);

   for (i = 0; i < nLibs; i++) {
      cur_bmaFile = bmas[i].bmaFile;
      sprintf_chk(personalBmaDataFile, MAX_FILE_PATH, "%s/%s_%d", base_dir, cur_bmaFile, MYTHREAD);
      get_rank_path(personalBmaDataFile, MYTHREAD);
#ifdef DEBUG
      fprintf(stderr, "Thread %d: Processing %s\n", MYTHREAD, personalBmaDataFile);
#endif
      cur_bma_fd = fopen_chk(personalBmaDataFile, "r");
      
      start_io = UPC_TICKS_NOW();
      
      fgets_result = fgets(line, MY_LINESIZE, cur_bma_fd);
      
      end_io = UPC_TICKS_NOW();
      io_read_time += UPC_TICKS_TO_SECS(end_io - start_io);
      int line_count = 0;
      
      
      while ( fgets_result  != NULL ) {
         line_count++;
         split_bma_entry(line, col0, col1, col2, col3, col4, col5);
         if ( (strcmp(col0, "SINGLE") == 0) && (strcmp(col1, "SPLINT") == 0) ) {
            processSplint(col2, col3, col4, col5, &bmas[i], local_index, local_buffs, link_heap);
         } else if (strcmp(col0, "PAIR") == 0) {
            processPair(col2, col3, col4, col5, &bmas[i], local_index, local_buffs, link_heap, &cpy_previousMaxInsertSize, minNetEndDistance);
         }
         
         start_io = UPC_TICKS_NOW();

         fgets_result = fgets(line, MY_LINESIZE, cur_bma_fd);
         
         end_io = UPC_TICKS_NOW();
         io_read_time += UPC_TICKS_TO_SECS(end_io - start_io);
      }
      
      fclose(cur_bma_fd);
#ifdef DEBUG
      fprintf(stderr, "Thread %d: Finished %s (%d lines)\n", MYTHREAD, personalBmaDataFile, line_count);
#endif
   }
   
   upc_barrier;
   serial_printf("Threads done processing links\n");
   
   /* Adds remaining links to the shared heaps */
   add_rest_links_to_shared_heaps(local_index, local_buffs, link_heap);
   
   /* Wait until all heaps have been fixed */
   upc_barrier;
   serial_printf("Threads done fixing heaps\n");
   
   endT = UPC_TICKS_NOW();
   double links_storing = UPC_TICKS_TO_SECS(endT - startT);
   
   serial_printf("Threads done processing the input files\n");

   startT = UPC_TICKS_NOW();

   
#ifdef DEBUG
   //UPC_ATOMIC_FADD_I64(&tot_links,link_heap.heap_indices[MYTHREAD] );
   //printf("Thread %d has %lld entries in its local stack\n", MYTHREAD,link_heap.heap_indices[MYTHREAD] );
   //upc_barrier;
   //if (MYTHREAD == 0) printf("Total links are: %lld \n", tot_links);
#endif
   
   /* Create local hash table by iterating on local heap and storing link_data to local_hashtable */
   datum_heap_pos = 0;
   link_heap_pos = 0;
   entries_in_my_stack = link_heap.heap_indices[MYTHREAD];
   create_local_heaps(entries_in_my_stack, &datum_heap, &linkList_heap);
   local_hashtable = create_link_hash_table(entries_in_my_stack);

   if(local_hashtable == NULL) 
       DIE("Could not create_link_hashtable %ld\n", entries_in_my_stack);
   assert(upc_threadof(link_heap.link_ptr + MYTHREAD) == MYTHREAD);
   my_stack = (link_t*) link_heap.link_ptr[MYTHREAD];
#ifdef DEBUG
   fprintf(stderr, "Thread %d: created link_hash_table: %ld\n", MYTHREAD, entries_in_my_stack);
#endif

   int64_t accepted_links = 0;
   int64_t nLinks = 0;
   FILE *myFD = NULL;
#ifdef DEBUG
   { 
       char myLogName[MAX_FILE_PATH];
       sprintf_chk(myLogName, MAX_FILE_PATH, "%s/linksLog_%d", base_dir, MYTHREAD);
       get_rank_path(myLogName, MYTHREAD);
       myFD = fopen_chk(myLogName, "w");
       if (myFD == NULL)
         DIE("Thread %d: Could not open the my logfile %s!\n", MYTHREAD, myLogName);
   } while (0);
#endif

   for (i=0; i < entries_in_my_stack; i++) {
      assert(datum_heap_pos < entries_in_my_stack);
      accepted_links += (int64_t) add_link_data(local_hashtable, (link_t*) (my_stack+i), &datum_heap_pos, datum_heap, &link_heap_pos, linkList_heap, bmas, &nLinks, myFD);
   }
   
#ifdef DEBUG
   fprintf(stderr, "Thread %d: added %ld links, accepted: %ld\n", MYTHREAD, nLinks, accepted_links);
#endif
   upc_barrier;
   serial_printf("Threads created local splints hash tables\n");
   //upc_barrier;
   
#ifdef DEBUG
   UPC_ATOMIC_FADD_I64(&tot_links,accepted_links);
   UPC_ATOMIC_FADD_I64(&tot_nLinks,nLinks);
   upc_barrier;
   if (MYTHREAD == 0) printf("Total accepted links are: %"PRId64" \n", tot_links);
   if (MYTHREAD == 0) printf("Total different link keys are: %"PRId64" \n", tot_nLinks);

#endif
   
   sprintf_chk(personalBmaDataFileOut, MAX_FILE_PATH, "%s/LINKS_OUTPUT_%d_%d", base_dir, linkUid, MYTHREAD);
   get_rank_path(personalBmaDataFileOut, MYTHREAD);
#ifdef DEBUG
   fprintf(stderr, "Thread %d: Writing links output: %s\n", MYTHREAD, personalBmaDataFileOut);
#endif
   outFD = fopen_chk(personalBmaDataFileOut, "w+");
   
   /*******************************************************/
   /* Now iterate over the hashtable and process the data */
   /*******************************************************/
   
   int64_t my_num_spans = 0;

   for (i=0; i < local_hashtable->size; i++) {
      //if (MYTHREAD == 0) printf("Index is %lld\n", i);
      cur_bucket = local_hashtable->table[i];
      while (cur_bucket != NULL) {
         
         minimumGap = - (merSize - 2);
         splintMaxDev = 2;
         spanMaxZ = 3;
         
         /* Set SPLINT frequency table to zero */
         memset(splints, 0, MAX_END_SEPARATION * sizeof(int64_t));
         /* Set span data structure to zero */
         for (j=0; j<nLibs; j++ ) {
            spans[j][0] = 0;
            spans[j][1] = 0;
            spans[j][2] = 0;
            spansGE[j] = 0.0;
         }
         nSplints = 0;
         nSpans = 0;
         
         /* Process the data in the the current link */
         cur_e1 = cur_bucket->end1_id;
         cur_e2 = cur_bucket->end2_id;
         cur_compressedTails = cur_bucket->compressedTails;
         if ( cur_compressedTails == 0 ) {
            tail1 =3;
            tail2 =3;
         }
         if ( cur_compressedTails == 1 ) {
            tail1 =3;
            tail2 =5;
         }
         if ( cur_compressedTails == 2 ) {
            tail1 =5;
            tail2 =3;
         }
         if ( cur_compressedTails == 3 ) {
            tail1 =5;
            tail2 =5;
         }
         
#ifdef DEBUG
         sprintf_chk(link, MY_LINESIZE, "Contig%d.%d<=>Contig%d.%d", cur_e1, tail1, cur_e2, tail2);
         fprintf(myFD, "%s\n", link);
         my_count = 0;
         my_checksum = 0;
         my_anomalies = 0;
#endif
         cur_type = cur_bucket->link_type;
         cur_data = cur_bucket->data;
         while (cur_data != NULL) {
            cur_lib = cur_data->lib_id;
            cur_endSep = cur_data->endSeparation;
            
#ifdef DEBUG
            my_count++;
            //fprintf(myFD, "%d, ", cur_endSep);
            my_checksum += cur_endSep;
#endif
            
            if ( cur_type == SPLINT ) {
               if ( cur_endSep >= (minimumGap - splintMaxDev) ) {
                  nSplints++;
                  splints[cur_endSep+OFFSET]++;
               }
            } else if ( cur_type == SPAN) {
               spanAnomaly = 0;
               if ( cur_endSep < minimumGap ) {
                  insertStdDev = bmas[cur_lib].stdDev;
                  spanZ = abs(minimumGap - cur_endSep)/ (1.0 * insertStdDev);
                  if (spanZ > spanMaxZ) {
#ifdef DEBUG
                     my_anomalies++;
#endif
                     spanAnomaly = 1;
                  }
               }
               if (spanAnomaly == 0) {
                  nSpans++;
                  spans[cur_lib][0] += 1;
                  spans[cur_lib][1] += cur_endSep;
               }
               
            }
            cur_data = cur_data->next;
         }
#ifdef DEBUG
         fprintf(myFD, "COUNT: %d, CHECKSUM: %d, ANOMALIES: %d, SPANS: %"PRId64"\n\n", my_count, my_checksum, my_anomalies, nSpans);
#endif
         /* END: Process the data in the the current link */
         
         /* Gap size as estimated from splints is taken to be the maximum frequency splint */
         splintGapEstimate = SPLINT_UNDEFINED;
         splintMaxFreq = 0;
         for (j=0; j<MAX_END_SEPARATION; j++) {
            if (splints[j] > splintMaxFreq ) {
               splintMaxFreq = splints[j];
               /* FIXME: Deal with negative values of endSeparation */
               splintGapEstimate = j - OFFSET;
            }
         }
         
         /* Gap size as estimated from spans is taken to be the weighted mean of spans */
         for (j=0; j < nLibs; j++) {
            nLibSpans = spans[j][0];
            if (nLibSpans > 0) {
               insertSize = bmas[j].insertSize;
               insertStdDev = bmas[j].stdDev;
               readLength = bmas[j].readLength;
               assert(readLength > 0);
               netSeparation = spans[j][1];
               meanGapEstimate = (1.0 * netSeparation/ (1.0 * nLibSpans));
               meanOffset = (double)insertSize - (double) meanGapEstimate;
               l1 = cur_bucket->o1_length;
               l2 = cur_bucket->o2_length;
#ifdef DEBUG
               sprintf_chk(link, MY_LINESIZE, "Contig%d.%d<=>Contig%d.%d", cur_e1, tail1, cur_e2, tail2);
               if (strcmp(link, "Contig1277.3<=>Contig5037.3") == 0) {
                  if (l1 < l2)
                     printf("Arguments are %f, %d, %d, %d, %d, %d, %d\n", meanOffset, merSize, readLength, l1, l2, insertSize, insertStdDev  );
                  else
                     printf("Arguments are %f, %d, %d, %d, %d, %d, %d\n", meanOffset, merSize, readLength, l2, l1, insertSize, insertStdDev  );
               }
#endif
               gapEstimate = (l1 < l2) ? estimateGapSize(meanOffset, merSize, readLength, l1, l2, 1.0*insertSize, 1.0*insertStdDev) : estimateGapSize(meanOffset, merSize, readLength, l2, l1, 1.0*insertSize, 1.0*insertStdDev) ;
               //spans[j][2] = gapEstimate;
               spansGE[j] = gapEstimate;
#ifdef DEBUG
               sprintf_chk(link, MY_LINESIZE, "Contig%d.%d<=>Contig%d.%d", cur_e1, tail1, cur_e2, tail2);
               if (strcmp(link, "Contig1277.3<=>Contig5037.3") == 0) {
                  printf("gap estimate returned is %f\n", gapEstimate  );
               }
#endif
            }
         }
         
         spanGapEstimate = 0.0;
         sumOfWeights = 0.0;
         spanGapUncertainty = -1;
         
         for (j=0; j < nLibs; j++) {
            if (!spans[j][0])
               continue;
            insertStdDev = bmas[j].stdDev;
            libWeight = (1.0*spans[j][0])/(1.0*insertStdDev*1.0*insertStdDev);
            sumOfWeights += 1.0*libWeight;
            //spanGapEstimate += (spans[j][2]) * libWeight;
            spanGapEstimate += (1.0 * spansGE[j]) * (1.0*libWeight);
         }
         
         if (sumOfWeights > EPSILON) {
            spanGapEstimate /= sumOfWeights;
            spanGapUncertaintyDouble = sqrt(1/sumOfWeights);
            spanGapUncertainty = 1;
         }
         
         anomalousSplints = 0;
         anomalousSpans = 0;
         
         /* Revisit data checking for consistency with estimate */
         cur_data = cur_bucket->data;
         while (cur_data != NULL) {
            cur_lib = cur_data->lib_id;
            cur_endSep = cur_data->endSeparation;
            if ( cur_type == SPLINT ) {
               if (splintGapEstimate != SPLINT_UNDEFINED) {
                  splintDev = abs(cur_endSep - splintGapEstimate);
                  if (splintDev > splintMaxDev) {
                     anomalousSplints++;
                  }
               }
            } else if ( cur_type == SPAN) {
               if (spanGapUncertainty != -1) {
                  insertStdDev = bmas[cur_lib].stdDev;
                  spanZ = fabs(cur_endSep - spanGapEstimate)/ (1.0*insertStdDev);
                  if (spanZ > spanMaxZ) {
#ifdef DEBUG
                     sprintf_chk(link, MY_LINESIZE, "Contig%d.%d<=>Contig%d.%d", cur_e1, tail1, cur_e2, tail2);
                     if (strcmp(link, "Contig1277.3<=>Contig5037.3") == 0) {
                        printf("SpanZ is %f , spanGapEstimate is %f\n", spanZ , spanGapEstimate );
                     }
#endif
                     anomalousSpans++;
                  }
               }
            }
            cur_data = cur_data->next;
         }
         
         if (splintGapEstimate != SPLINT_UNDEFINED) {
            //strcpy(link, "Contig");
            //myitoa(cur_e1, temp_string, 10);
            //strcpy(link, temp_string);
            //strcat(link, "<=>");
            //strcpy(link, "Contig");
            //myitoa(cur_e2, temp_string, 10);
            //strcpy(link, temp_string);
            myPRINTlinks++;
            sprintf_chk(link, MY_LINESIZE, "Contig%d.%d<=>Contig%d.%d", cur_e1, tail1, cur_e2, tail2);
            fprintf(outFD, "SPLINT\t%s\t%ld|%ld|%ld\t%d\n", link, splintMaxFreq, anomalousSplints, nSplints, splintGapEstimate);
         }
         
         if (spanGapUncertainty != -1) {
            //strcpy(link, "Scaffold");
            //myitoa(cur_e1, temp_string, 10);
            //strcpy(link, temp_string);
            //strcat(link, "<=>");
            //strcpy(link, "Scaffold");
            //myitoa(cur_e2, temp_string, 10);
            //strcpy(link, temp_string);
            myPRINTlinks++;
            if (scaffoldSpans == 0) {
               sprintf_chk(link, MY_LINESIZE, "Contig%d.%d<=>Contig%d.%d", cur_e1, tail1, cur_e2, tail2);
            } else {
               sprintf_chk(link, MY_LINESIZE, "Scaffold%d.%d<=>Scaffold%d.%d", cur_e1, tail1, cur_e2, tail2);
            }
            fprintf(outFD, "SPAN\t%s\t%ld|%ld\t%.0f\t%.0f\n", link, anomalousSpans, nSpans, spanGapEstimate, spanGapUncertaintyDouble);
            my_num_spans++;
         }
         
         cur_bucket = cur_bucket->next;
      }
   }


#ifdef DEBUG
   fclose(myFD);
#endif
   /* TODO:  Add functionality to add result at separate files  */
   
   upc_barrier;
   end = UPC_TICKS_NOW();
   
   endT = UPC_TICKS_NOW();
   double processing_links = UPC_TICKS_TO_SECS(endT - startT);
   serial_printf("Threads done with output files\n");
   
   if (MYTHREAD == 0) {
      printf("\nTime for executing bmaToLinks : %d seconds\n" ,((int)UPC_TICKS_TO_SECS(end-start)));
      printf("Time for storing links data (includes communication) is %f seconds\n", links_storing);
      printf("Time for processing links is %f seconds\n", processing_links);
      printf("I/O read time is %f seconds\n", io_read_time);


   }

   bupc_atomicI64_fetchadd_strict(&PRINTnLinks, myPRINTlinks);
   upc_barrier;
   
   /* Report summary statistics for each library */
   if ( MYTHREAD == 0 ) {
      printf("\n************* Report summary statistics for each library ***********\n");
   }
   for (i = 0; i < nLibs; i++) {
      upc_barrier;
      nTRUNC = 0;
      nSMALL = 0;
      nINORM = 0;
      nINNIE = 0;
      nSELFL = 0;
      nEDIST = 0;
      nREDUN = 0;
      nACCPT = 0;
      nISHIRT = 0;
      upc_barrier;
      
      bupc_atomicI64_fetchadd_strict(&nTRUNC, bmas[i].libPairSummary[TRUNCA]);
      bupc_atomicI64_fetchadd_strict(&nINNIE, bmas[i].libPairSummary[INNIE]);
      bupc_atomicI64_fetchadd_strict(&nINORM, bmas[i].libPairSummary[INORM]);
      bupc_atomicI64_fetchadd_strict(&nREDUN, bmas[i].libPairSummary[REDUN]);
      bupc_atomicI64_fetchadd_strict(&nSELFL, bmas[i].libPairSummary[SELFL]);
      bupc_atomicI64_fetchadd_strict(&nSMALL, bmas[i].libPairSummary[SMALL]);
      bupc_atomicI64_fetchadd_strict(&nISHIRT, bmas[i].libPairSummary[ISHIRT]);
      bupc_atomicI64_fetchadd_strict(&nACCPT, bmas[i].libPairSummary[ACCPT]);
      bupc_atomicI64_fetchadd_strict(&nEDIST, bmas[i].libPairSummary[EDIST]);
      bupc_atomicI64_fetchadd_strict(&tot_spans, my_num_spans / nLibs);
      
      upc_barrier;

      if ( MYTHREAD == 0 ) {
         total = nTRUNC + nSMALL + nINORM + nINNIE + nSELFL + nEDIST + nREDUN + nACCPT + nISHIRT;
         if (total != 0) {
            printf("LIB: %s\n", bmas[i].bmaFile);
            printf("Number of spans: %ld\n", tot_spans);
            printf("Percentage of TRUNC: %.3f %% \n", (1.0 * nTRUNC) / (1.0 * total ) * 100.0 );
            printf("Percentage of SMALL: %.3f %% \n", (1.0 * nSMALL) / (1.0 * total ) * 100.0 );
            printf("Percentage of INORM: %.3f %% \n", (1.0 * nINORM) / (1.0 * total ) * 100.0 );
            printf("Percentage of ISHIRT: %.3f %% \n", (1.0 * nISHIRT) / (1.0 * total ) * 100.0 );
            printf("Percentage of INNIE: %.3f %% \n", (1.0 * nINNIE) / (1.0 * total ) * 100.0 );
            printf("Percentage of SELFL: %.3f %% \n", (1.0 * nSELFL) / (1.0 * total ) * 100.0 );
            printf("Percentage of EDIST: %.3f %% \n", (1.0 * nEDIST) / (1.0 * total ) * 100.0 );
            printf("Percentage of REDUN: %.3f %% \n", (1.0 * nREDUN) / (1.0 * total ) * 100.0 );
            printf("Percentage of ACCPT: %.3f %% \n", (1.0 * nACCPT) / (1.0 * total ) * 100.0 );
            printf("Total links processed: %ld\n", total);
            printf("********************************************\n\n");
         }
      }
      upc_barrier;
   }
   
   if (MYTHREAD == 0) {
      char countFileName[MAX_FILE_PATH];
      sprintf_chk(countFileName, MAX_FILE_PATH, "nLinks-%d.txt", linkUid);
      get_rank_path(countFileName, -1);
      FILE *countFD = fopen_chk(countFileName, "w+" );

      fprintf(countFD, "%ld\n", PRINTnLinks);
      fclose(countFD);
      
      sprintf_chk(countFileName, MY_LINESIZE, "linksMeta-%d", linkUid);
      get_rank_path(countFileName, -1);
      countFD = fopen_chk(countFileName, "w+" );
      /* use that last library's insertSize & stdDev */
      fprintf(countFD, "%d\t%d\tLINKS_OUTPUT_%d\n", bmas[nLibs-1].insertSize, bmas[nLibs-1].stdDev, linkUid);
      fclose(countFD);
   }
   
   upc_barrier;

	if (!MYTHREAD)
		printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
			   ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
   
   return 0;
}
