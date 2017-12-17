#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <ctype.h>
#include <upc.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <libgen.h>

#include "../../common/optlist.h"
#include "../../common/common.h"
#include "../../common/upc_compatibility.h"
#include "analyzerUtils.h"
#include "locationHash.h"

shared int64_t nLoci = 0;
shared int64_t redundancy = 0;

shared int64_t nFullPairs = 0;
shared int64_t nNonStandard = 0;
int64_t myNonStandard = 0;
int64_t meanRedundancy = 0;

shared int64_t num_found_fails = 0;
shared int64_t num_length_fails = 0;
shared int64_t num_tests = 0;

int main(int argc, char **argv) {
   upc_tick_t start_time = upc_ticks_now();
   int64_t myFullPairs = 0;
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "l:i:s:m:M:RAE:L:F:T:U:t:c:B:N:Z:G:C:p:");
   int reverseComplement = 0, innieRemoval = 0, endAversion = 0 , minTestLength = 0, minFreqReported = 10, fivePrimeWiggleRoom = 5, threePrimeWiggleRoom = 5, truncate = 0, minMatch = 0, estimatedInsertSize = 0, estimatedSigma = 0, testMode = 1, insertSize = 0, insertSigma = 0, coresPerNode = 1;
   char *libname = NULL;
   char outputMerAligner[MAX_FILE_PATH];
   FILE *laneFD1, *laneFD2;
   int64_t i, j, k, pos1, pos2, found_flag, partial_res, test1, test2, test3, strand1, strand2, sStart1, sStart2, sStop1, sStop2, sLength1, locs[4], min, max, extent, safeDistance, leftBoundary, rightBoundary, orientation, qStart1, qStart2, qStop1, qStop2, qLength2, qLength1, sLength2;
   align_info result1, result2;
   char *lineBuffers = (char*) calloc_chk(MAX_LINE_SIZE * 10, 1);
   char *alignmentBuffer1 = lineBuffers;
   char *copyAlignmentBuffer1 = alignmentBuffer1 + MAX_LINE_SIZE;
   char *firstAlignmentBuffer1 = copyAlignmentBuffer1 + MAX_LINE_SIZE;
   char *alignmentBuffer2 = firstAlignmentBuffer1 + MAX_LINE_SIZE;
   char *copyAlignmentBuffer2 = alignmentBuffer2 + MAX_LINE_SIZE;
   char *firstAlignmentBuffer2 = copyAlignmentBuffer2 + MAX_LINE_SIZE;
   char *readIdLane1 = firstAlignmentBuffer2 + MAX_LINE_SIZE;
   char *newReadIdLane1 = readIdLane1 + MAX_LINE_SIZE;
   char *readIdLane2 = newReadIdLane1 + MAX_LINE_SIZE;
   char *newReadIdLane2 = readIdLane2 + MAX_LINE_SIZE;
   char *resRead1, *resRead2;
   int validAlign1, validAlign2;
   UPC_TICK_T start, end;
   int coresToUsePerNode = 1;
   int shortPair = 600, trimOff;
   int threshold = 10, subjectID;
   int count;
   int64_t totalAlignments = 10;
   hash_table_t *locationHashTable;
   memory_heap_t memory_heap;
   bucket_t *cur_bucket;
   location_t *cur_location;
   int64_t hash_table_size, myLoci = 0, myRedundancy = 0, mean, meanSquare, n, size, freq;
   double meanDouble, stdDev;
   const char *base_dir = ".";
   
   int srfInput = 0;
   char *srfSuffixName;
   int totalScaffolds, totalContigs;
   
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'Z':
            /* FIXME: Standardize srf file output name conventions */
            srfInput = 1;
            srfSuffixName = thisOpt->argument;
            break;
         case 'G':
            totalScaffolds = atoi(thisOpt->argument);
            break;
         case 'C':
            totalContigs = atoi(thisOpt->argument);
            break;
         case 'R':
            reverseComplement = 1;
            break;
         case 'A':
            innieRemoval = 1;
            break;
         case 'l':
            libname = thisOpt->argument;
            break;
         case 'E':
            endAversion = atoi(thisOpt->argument);
            break;
         case 'L':
            minTestLength = atoi(thisOpt->argument);
            break;
         case 'U':
            truncate = atoi(thisOpt->argument);
            break;
         case 'T':
            threePrimeWiggleRoom = atoi(thisOpt->argument);
            break;
         case 't':
            threshold = atoi(thisOpt->argument);
            break;
         case 'M':
            minFreqReported = atoi(thisOpt->argument);
            break;
         case 'F':
            fivePrimeWiggleRoom = atoi(thisOpt->argument);
            break;
         case 'm':
            minMatch = atoi(thisOpt->argument);
            break;
         case 'i':
            estimatedInsertSize = atoi(thisOpt->argument);
            break;
         case 's':
            estimatedSigma = atoi(thisOpt->argument);
            break;
         case 'p':
            coresToUsePerNode = atoi(thisOpt->argument);
            break;
         case 'c':
            totalAlignments = __atoi64(thisOpt->argument);
         case 'B':
            base_dir = thisOpt->argument;
            break;
         case 'N':
            coresPerNode = atoi(thisOpt->argument);
            break;
         default:
            break;
      }
      free(thisOpt);
   }
   
   if (libname == NULL)
      DIE("You must specify a -l libname\n");
   
   start = UPC_TICKS_NOW();
   
   int64_t my_num_tests = 0, my_num_length_fails = 0, my_num_found_fails = 0, my_num_valid = 0;
   int64_t *local_histo = (int64_t*) calloc_chk(HISTOGRAM_SIZE, sizeof(int64_t));
   
   /* Create Location hash table and relevant memory heap */
   upc_barrier;
   locationHashTable = create_hash_table(totalAlignments + 200000*THREADS, &memory_heap);
   upc_barrier;
   
   shared[1] int64_t *scaffLengths;
   shared[1] contigScaffoldMap_t *contigScaffoldMap;
   contigScaffoldMap_t contigScaffoldMapEntry;
   char srfLine[MAX_LINE_SIZE];
   char srfFileName[MAX_FILE_PATH];
   FILE *srfFD;
   char *fgets_result;
   int splitRes, scaffoldId, contigId, cStrand, sStart, sEnd, subject1, subject2;
   upc_tick_t atomics_timer;
   
   /******************************/
   /*  Read scaffold report file */
   /******************************/
   
   if (srfInput) {
      if (!srfSuffixName)
         SDIE("You must specify a -Z srfSuffixName\n");
      sprintf(srfFileName, "%s/%s_%d", base_dir, srfSuffixName, MYTHREAD);
      srfFD = fopen_rank_path(srfFileName, "r", MYTHREAD);
      if (!totalScaffolds || !totalContigs)
         SDIE("You must specify a -C totalContigs and a -G totalScaffolds\n");
      
      /* Build scaffLengths data structure */
      scaffLengths = (shared[1] int64_t*) upc_all_alloc(totalScaffolds+1, sizeof(int64_t));
      if (scaffLengths == NULL)
         DIE("Could not allocate %lld scaffolds!\n", (long long) totalScaffolds);
      for (i = MYTHREAD; i < totalScaffolds; i += THREADS) {
         scaffLengths[i] = 0;
      }
      
      /* Build contigScaffoldMap data structure */
      contigScaffoldMap = (shared[1] contigScaffoldMap_t*) upc_all_alloc(totalContigs+1,sizeof(contigScaffoldMap_t));
      if (contigScaffoldMap == NULL)
         DIE("Could not allocate %lld totalContigs!\n", (long long) totalContigs);
      for (i = MYTHREAD; i < totalContigs; i += THREADS) {
         contigScaffoldMap[i].scaffID = UNDEFINED;
      }
      
      upc_barrier;
      
      fgets_result = fgets(srfLine, MAX_LINE_SIZE, srfFD);
      
      while ( fgets_result  != NULL) {
         splitRes = splitSrfLine(srfLine, &scaffoldId, &contigId, &sStart, &sEnd, &cStrand);
         if (splitRes == SUCCESS) {
            
            /* Update scaffsLength array */
            if (scaffLengths[scaffoldId] == 0) {
               scaffLengths[scaffoldId] = sEnd;
            } else {
               if (sEnd > scaffLengths[scaffoldId]) {
                  scaffLengths[scaffoldId] = sEnd;
               }
            }
            
            /* Update contigScaffoldMap data structure */
            contigScaffoldMapEntry.cStrand = cStrand;
            contigScaffoldMapEntry.scaffID = scaffoldId;
            contigScaffoldMapEntry.sStart = sStart;
            contigScaffoldMapEntry.sEnd = sEnd;
            contigScaffoldMap[contigId] = contigScaffoldMapEntry;
         }
         
         fgets_result = fgets(srfLine, MAX_LINE_SIZE, srfFD);
      }
      
      fclose(srfFD);
   }
   
   upc_barrier;
   
   /* Allocation of global histogram */
   shared[1] int64_t *shared_histo = (shared[1] int64_t*) upc_all_alloc(HISTOGRAM_SIZE, sizeof(int64_t));
   
   upc_forall (i=0; i<HISTOGRAM_SIZE; i++; i) {
      shared_histo[i] = 0;
   }
   
   upc_barrier;
   
   int numFilesPerThread = coresPerNode / coresToUsePerNode;

   if (MYTHREAD % numFilesPerThread == 0) {
       //printf("Thread %d is computing with %d files\n", MYTHREAD, numFilesPerThread);
      
      int scaffLength1, scaffLength2, contigScaffStrand1, contigScaffStrand2, scaffStrand1, 
         scaffStrand2, scaffStart1, scaffStart2, scaffStop1, scaffStop2, sLengthCutoff;
      
      for (int fileIndex = MYTHREAD; fileIndex < MYTHREAD + numFilesPerThread; fileIndex++) {
         sprintf(outputMerAligner, "%s/%s-merAlignerOutput_%d_Read1", base_dir, libname, fileIndex);
         laneFD1 = fopen_rank_path(outputMerAligner, "r", fileIndex);
         sprintf(outputMerAligner, "%s/%s-merAlignerOutput_%d_Read2", base_dir, libname, fileIndex);
         laneFD2 = fopen_rank_path(outputMerAligner, "r", fileIndex);
         //printf("Thread %d is reading %s\n", MYTHREAD, outputMerAligner);
   
         resRead1 = fgets(firstAlignmentBuffer1, MAX_LINE_SIZE, laneFD1);
         resRead2 = fgets(firstAlignmentBuffer2, MAX_LINE_SIZE, laneFD2);
   
         while (myFullPairs < totalAlignments) {
      
            if (resRead1 == NULL && resRead2 == NULL) break;
      
            /* Parse the first alignment of read (lane 1) */
            if (resRead1 != NULL) {
               validAlign1 = splitAlignment(firstAlignmentBuffer1, readIdLane1, &(result1.aligned_read_locs[0]), &(result1.aligned_read_locs[1]), &(result1.read_length), &(result1.aligned_contig_ids[0]), &(result1.aligned_contig_locs[0]), &(result1.aligned_contig_locs[1]), &(result1.aligned_contig_lengths[0]), &(result1.aligned_strand[0]));
               if (validAlign1 == 0) {
                  result1.contigs_matched = 0;
               } else {
                  result1.contigs_matched = 1;
                  my_num_valid++;
               }
            }
      
            /* Read the the rest alignments of read (lane 1) */
            if (resRead1 != NULL) {
               resRead1 = fgets(alignmentBuffer1, MAX_LINE_SIZE, laneFD1);
               memcpy(copyAlignmentBuffer1, alignmentBuffer1, MAX_LINE_SIZE * sizeof(char));
               if (resRead1 != NULL) {
                  validAlign1 = splitAlignment(copyAlignmentBuffer1, newReadIdLane1, &(result1.aligned_read_locs[2*result1.contigs_matched]), &(result1.aligned_read_locs[2*result1.contigs_matched+1]), &(result1.read_length), &(result1.aligned_contig_ids[result1.contigs_matched]), &(result1.aligned_contig_locs[2*result1.contigs_matched]), &(result1.aligned_contig_locs[2*result1.contigs_matched+1]), &(result1.aligned_contig_lengths[result1.contigs_matched]), &(result1.aligned_strand[result1.contigs_matched]));
            
                  while ( strcmp(readIdLane1, newReadIdLane1) == 0 ) {
                     result1.contigs_matched++;
                     resRead1 = fgets(alignmentBuffer1, MAX_LINE_SIZE, laneFD1);
                     memcpy(copyAlignmentBuffer1, alignmentBuffer1, MAX_LINE_SIZE * sizeof(char));
                     if (resRead1 == NULL) break;
                     validAlign1 = splitAlignment(copyAlignmentBuffer1, newReadIdLane1, &(result1.aligned_read_locs[2*result1.contigs_matched]), &(result1.aligned_read_locs[2*result1.contigs_matched+1]), &(result1.read_length), &(result1.aligned_contig_ids[result1.contigs_matched]), &(result1.aligned_contig_locs[2*result1.contigs_matched]), &(result1.aligned_contig_locs[2*result1.contigs_matched+1]), &(result1.aligned_contig_lengths[result1.contigs_matched]), &(result1.aligned_strand[result1.contigs_matched]));
                  }
            
                  if (resRead1 != NULL) {
                     memcpy(firstAlignmentBuffer1, alignmentBuffer1, MAX_LINE_SIZE * sizeof(char) );
                     strcpy(readIdLane1, newReadIdLane1);
                  }
               }
            }
      
            /* Parse the first alignment of read (lane 2) */
      
            if (resRead2 != NULL) {
               validAlign2 = splitAlignment(firstAlignmentBuffer2, readIdLane2, &(result2.aligned_read_locs[0]), &(result2.aligned_read_locs[1]), &(result2.read_length), &(result2.aligned_contig_ids[0]), &(result2.aligned_contig_locs[0]), &(result2.aligned_contig_locs[1]), &(result2.aligned_contig_lengths[0]), &(result2.aligned_strand[0]));
               if (validAlign2 == 0) {
                  result2.contigs_matched = 0;
               } else {
                  result2.contigs_matched = 1;
               }
            }
      
            /* Read the the rest alignments of read (lane 2) */
            if (resRead2 != NULL) {
               resRead2 = fgets(alignmentBuffer2, MAX_LINE_SIZE, laneFD2);
               memcpy(copyAlignmentBuffer2, alignmentBuffer2, MAX_LINE_SIZE * sizeof(char));
         
               if (resRead2 != NULL) {
                  validAlign2 = splitAlignment(copyAlignmentBuffer2, newReadIdLane2, &(result2.aligned_read_locs[2*result2.contigs_matched]), &(result2.aligned_read_locs[2*result2.contigs_matched+1]), &(result2.read_length), &(result2.aligned_contig_ids[result2.contigs_matched]), &(result2.aligned_contig_locs[2*result2.contigs_matched]), &(result2.aligned_contig_locs[2*result2.contigs_matched+1]), &(result2.aligned_contig_lengths[result2.contigs_matched]), &(result2.aligned_strand[result2.contigs_matched]));
            
                  while ( strcmp(readIdLane2, newReadIdLane2) == 0 ) {
                     result2.contigs_matched++;
                     resRead2 = fgets(alignmentBuffer2, MAX_LINE_SIZE, laneFD2);
                     memcpy(copyAlignmentBuffer2, alignmentBuffer2, MAX_LINE_SIZE * sizeof(char));
                     if (resRead2 == NULL) break;
                     validAlign2 = splitAlignment(copyAlignmentBuffer2, newReadIdLane2, &(result2.aligned_read_locs[2*result2.contigs_matched]), &(result2.aligned_read_locs[2*result2.contigs_matched+1]), &(result2.read_length), &(result2.aligned_contig_ids[result2.contigs_matched]), &(result2.aligned_contig_locs[2*result2.contigs_matched]), &(result2.aligned_contig_locs[2*result2.contigs_matched+1]), &(result2.aligned_contig_lengths[result2.contigs_matched]), &(result2.aligned_strand[result2.contigs_matched]));
                  }
            
                  if (resRead2 != NULL) {
                     memcpy(firstAlignmentBuffer2, alignmentBuffer2, MAX_LINE_SIZE * sizeof(char) );
                     strcpy(readIdLane2, newReadIdLane2);
                  }
               }
            }
      
            /*****************************/
            /* Find insert size estimate */
            /*****************************/
      
            if ( (result1.contigs_matched != 1) || (result2.contigs_matched != 1) ) {
               continue;
            }
      
            my_num_tests++;
      
            /* Extract parameters of the alignments */
            strand1 = result1.aligned_strand[0];
            sStart1 = result1.aligned_contig_locs[0];
            qStart1 = result1.aligned_read_locs[0];
            sStop1 = result1.aligned_contig_locs[1];
            qStop1 = result1.aligned_read_locs[1];
            qLength1 = result1.read_length;
            sLength1 = result1.aligned_contig_lengths[0];
            strand2 = result2.aligned_strand[0];
            sStart2 = result2.aligned_contig_locs[0];
            qStart2 = result2.aligned_read_locs[0];
            sStop2 = result2.aligned_contig_locs[1];
            qStop2 = result2.aligned_read_locs[1];
            sLength2 = result2.aligned_contig_lengths[0];
            qLength2 = result2.read_length;
            subject1 = result1.aligned_contig_ids[0];
            subject2 = result2.aligned_contig_ids[0];

            /* Apply truncation conditions to lane 1 read */
            if (truncate == 5) {
               if ( (qStart1 > 1) || (qStop1 < minMatch) ) {
                  continue;
               }
               trimOff = qStop1 - minMatch;
               qStop1 = minMatch;
               qLength1 = minMatch;
               if (strand1 == PLUS) {
                  sStop1 -= trimOff;
               } else {
                  sStart1 += trimOff;
               }
            } else if (truncate == 3) {
               if ( (qLength1 - qStart1 + 1 < minMatch) || (qStop1 < qLength1) ) {
                  continue;
               }
               trimOff = (qLength1 - minMatch + 1) - qStart1;
               qStart1 = 1;
               qStop1 = minMatch;
               qLength1 = minMatch;
               if (strand1 == PLUS) {
                  sStart1 += trimOff;
               } else {
                  sStop1 -= trimOff;
               }
            }
      
            /* Apply truncation conditions to lane 2 read */
            if (truncate == 5) {
               if ( (qStart2 > 1) || (qStop2 < minMatch) ) {
                  continue;
               }
               trimOff = qStop2 - minMatch;
               qStop2 = minMatch;
               qLength2 = minMatch;
               if (strand2 == PLUS) {
                  sStop2 -= trimOff;
               } else {
                  sStart2 += trimOff;
               }
            } else if (truncate == 3) {
               if ( (qLength2 - qStart2 + 1 < minMatch) || (qStop2 < qLength2) ) {
                  continue;
               }
               trimOff = (qLength2 - minMatch + 1) - qStart2;
               qStart2 = 1;
               qStop2 = minMatch;
               qLength2 = minMatch;
               if (strand2 == PLUS) {
                  sStart2 += trimOff;
               } else {
                  sStop2 -= trimOff;
               }
            }
      
            /* Translate locations if there is an srf input file */
            if (srfInput) {
               contigScaffoldMapEntry = contigScaffoldMap[subject1];
               if (contigScaffoldMapEntry.scaffID == UNDEFINED) {
                  continue;
               }
         
               scaffLength1 = scaffLengths[contigScaffoldMapEntry.scaffID];
               contigScaffStrand1 = contigScaffoldMapEntry.cStrand;
               scaffStrand1 = (contigScaffStrand1 == strand1) ? PLUS : MINUS;
               scaffStart1 = (contigScaffStrand1 == PLUS) ? contigScaffoldMapEntry.sStart + sStart1 -1 : contigScaffoldMapEntry.sEnd - sStop1 +1 ;
               scaffStop1 = (contigScaffStrand1 == PLUS) ? contigScaffoldMapEntry.sStart + sStop1 -1 : contigScaffoldMapEntry.sEnd - sStart1 +1 ;
         
               subject1 = contigScaffoldMapEntry.scaffID;
               sLength1 = scaffLength1;
               strand1 = scaffStrand1;
               sStart1 = scaffStart1;
               sStop1 = scaffStop1;
         
               contigScaffoldMapEntry = contigScaffoldMap[subject2];
               if (contigScaffoldMapEntry.scaffID == UNDEFINED) {
                  continue;
               }
         
               scaffLength2 = scaffLengths[contigScaffoldMapEntry.scaffID];
               contigScaffStrand2 = contigScaffoldMapEntry.cStrand;
               scaffStrand2 = (contigScaffStrand2 == strand2) ? PLUS : MINUS;
               scaffStart2 = (contigScaffStrand2 == PLUS) ? contigScaffoldMapEntry.sStart + sStart2 -1 : contigScaffoldMapEntry.sEnd - sStop2 +1 ;
               scaffStop2 = (contigScaffStrand2 == PLUS) ? contigScaffoldMapEntry.sStart + sStop2 -1 : contigScaffoldMapEntry.sEnd - sStart2 +1 ;
         
               subject2 = contigScaffoldMapEntry.scaffID;
               sLength2 = scaffLength2;
               strand2 = scaffStrand2;
               sStart2 = scaffStart2;
               sStop2 = scaffStop2;
            }
      
            /* If the read does not align to the same pair, then reject it */
            if ( subject1 != subject2 ) {
               continue;
            }
      
            /* endAversion test */
            if (endAversion) {
               if ( (sStop1 <= endAversion) || (sStart1 >= sLength1 - endAversion) ) {
                  continue;
               }
               if ( (sStop2 <= endAversion) || (sStart2 >= sLength2 - endAversion) ) {
                  continue;
               }
            }
      
            if (testMode) {
               sLengthCutoff = 2 * estimatedInsertSize + 8 * estimatedSigma;
               if (sLength1 <= sLengthCutoff) {
                  my_num_length_fails++;
                  continue;
               }
            }
      
            if (reverseComplement) {
               strand1 = (strand1 == PLUS) ? MINUS : PLUS;
               strand2 = (strand2 == PLUS) ? MINUS : PLUS;
            }
      
            test1 = testAlignStatus(qStart1, qStop1, qLength1, sStart1, sStop1, sLength1, strand1, fivePrimeWiggleRoom, threePrimeWiggleRoom, innieRemoval, shortPair);
      
            if (test1 == FAIL) continue;
      
            test2 = testAlignStatus(qStart2, qStop2, qLength2, sStart2, sStop2, sLength2, strand2, fivePrimeWiggleRoom, threePrimeWiggleRoom, innieRemoval, shortPair);
      
            if (test2 == FAIL) continue;
      
      
            if (strand1 == strand2) {
               test3 = FAIL;
               myNonStandard++;
            } else {
               test3 = SUCCESS;
            }
      
            /* Test orientation */
      
            orientation = UNDEFINED;
            if ( (test1 == SUCCESS) && (test2 == SUCCESS) && (test3 == SUCCESS) ) {
               if ( ( (strand1 == PLUS) && (sStart1 < sStart2) && (sStart1 < sStop2) ) || ( (strand2 == PLUS) && (sStart2 < sStart1) && (sStart2 < sStop1) ) ) {
                  orientation = CONVERGENT;
               } else if ( ( (strand1 == PLUS) && (sStart1 > sStop2) ) || ( (strand2 == PLUS) && (sStart2 > sStop1) ) ) {
                  orientation = DIVERGENT;
               }
         
               if (orientation != CONVERGENT) {
                  myNonStandard++;
                  test3 = FAIL;
               }
            }
      
            /* Test boundaries */
            if ( (test1 == SUCCESS) && (test2 == SUCCESS) && (test3 == SUCCESS) ) {
               locs[0] = sStart1;
               locs[1] = sStart2;
               locs[2] = sStop1;
               locs[3] = sStop2;
               qsort(locs, 4, sizeof(int64_t), cmpfunc);
               min = locs[0];
               max = locs[3];
               safeDistance = estimatedInsertSize+4*estimatedSigma;
               leftBoundary = safeDistance;
               rightBoundary = sLength1 - safeDistance;
               if ( (max < leftBoundary) || (min > rightBoundary) ) {
                  test3 = FAIL;
               }
            }
      
            /* Test location */
            if ( (test1 == SUCCESS) && (test2 == SUCCESS) && (test3 == SUCCESS) ) {
               subjectID = subject1;
               count = lookupAndIncrementLocation(subjectID, min, max, locationHashTable, &memory_heap);
               if (count != 0) {
                  test1 = FAIL;
               }
            }
      
            /* Test innie removal */
            if ( (test1 == SUCCESS) && (test2 == SUCCESS) && (test3 == SUCCESS) ) {
               extent = max-min+1;
               if (innieRemoval) {
                  if (extent < shortPair) {
                     myNonStandard++;
                     test1 = FAIL;
                  }
               }
            }
      
            /* Have passed all the tests, so add this entry in the histogram */
            if ( (test1 == SUCCESS) && (test2 == SUCCESS) && (test3 == SUCCESS) ) {
               extent = max-min+1;
               myFullPairs++;
               /* Check for bounds in histogram */
               if ( extent < HISTOGRAM_SIZE )
                  local_histo[extent]++;
            }
      
         } // while ( moreAlignmentsExist ) { 
      
         fclose(laneFD1);
         fclose(laneFD2);

         if (myFullPairs == totalAlignments) break;
      } // file loop
      //if (myFullPairs == totalAlignments) 
      //   printf("Thread %d: Found an adequate sample: %ld\n", MYTHREAD, myFullPairs);
      //else 
      //   printf("WARNING: Thread %d did not find sufficient data ( %ld < %ld ) in files for a full sample\n", 
      //          MYTHREAD, myFullPairs, totalAlignments);

   } //   if ((MYTHREAD % coresPerNode) < coresToUsePerNode) {
   
   atomics_timer = UPC_TICKS_NOW();
   for (i=0; i<HISTOGRAM_SIZE; i++) {
      if (local_histo[i] !=0 ) {
         partial_res = bupc_atomicI64_fetchadd_strict(&shared_histo[i], local_histo[i]);
      }
   }
      
   upc_barrier;
   
   if (MYTHREAD == 0) {
      for (i = 0; i < HISTOGRAM_SIZE; i++) {
         local_histo[i] = shared_histo[i];
      }
      end = UPC_TICKS_NOW();
      printf("\nTime for computing histogrammer : %.2f s, reduction %.2f s\n",
             (UPC_TICKS_TO_SECS(end-start)), (UPC_TICKS_TO_SECS(end-atomics_timer)));
   }
   
   upc_barrier;
   
   /* Report results in a more formal way */
   /* Each thread iterates overs its entries of testLocations table to measure nLoci and redundancy */
   hash_table_size = locationHashTable->size;
   for (i=MYTHREAD; i < hash_table_size; i+=THREADS) {
      cur_location = (location_t*) locationHashTable->table[i].head;
      while (cur_location != NULL) {
         myLoci++;
         myRedundancy += cur_location->count;
         cur_location = (location_t*) cur_location->next;
      }
   }
   
   
   bupc_atomicI64_fetchadd_strict(&nLoci, myLoci);
   bupc_atomicI64_fetchadd_strict(&redundancy, myRedundancy);
   bupc_atomicI64_fetchadd_strict(&nFullPairs, myFullPairs);
   bupc_atomicI64_fetchadd_strict(&nNonStandard, myNonStandard);
   
   bupc_atomicI64_fetchadd_strict(&num_tests, my_num_tests);
   bupc_atomicI64_fetchadd_strict(&num_length_fails, my_num_length_fails);
   
   upc_barrier;
   
   /* Thread 0 will report the statistics */
   
   if (MYTHREAD == 0) {
      printf("# Carried out %ld tests; failures: too short %ld (%.3f)\n",
             num_tests, num_length_fails, (double)num_length_fails / num_tests);
      printf("# %ld sampled canonical pairs (%ld non-standard pairs discarded); mean redundancy = %.1f\n", nFullPairs, nNonStandard ,(1.0 * redundancy) / (1.0 * nLoci));
      
      mean = 0;
      meanSquare = 0;
      n = 0;
      for (i = 0; i < HISTOGRAM_SIZE; i++) {
         size = i;
         if (size >= minTestLength) {
            freq = local_histo[size];
            if (freq >= minFreqReported) {
               mean += freq * size;
               meanSquare += freq * size * size;
               n += freq;
            }
         }
      }
      
      meanDouble = (1.0 * mean)/(1.0 * n);
      stdDev = sqrt((1.0 * meanSquare)/(1.0 * n) - meanDouble * meanDouble);
      if(isnan(meanDouble) || isnan(stdDev)) {
         printf("Warning: Could not determine mean insert size, using estimated (mean=%ld / n=%ld)\n", mean, n);
         meanDouble = estimatedInsertSize;
         stdDev = estimatedSigma;
      }
      
/**** CHECKING */
//      meanDouble = estimatedInsertSize;
//      stdDev = estimatedSigma;
/****/
      printf("# %.0f +/- %.0f (if this is very different from %d +/- %d consider rerun!)\n",
             meanDouble, stdDev, estimatedInsertSize, estimatedSigma);
      printf("# Pair separation distribution:\n");
      
      
      char countFileName[MAX_FILE_PATH];
      sprintf(countFileName,"%s-insert.txt", libname);
      get_rank_path(countFileName, -1);
      FILE *countFD = fopen_chk(countFileName, "w+" );
      fprintf(countFD, "%.0f\n", meanDouble);
      fclose(countFD);
      
      sprintf(countFileName,"%s-std.txt", libname);
      get_rank_path(countFileName, -1);
      countFD = fopen_chk(countFileName, "w+" );
      fprintf(countFD, "%.0f\n", stdDev);
      fclose(countFD);
   }

   
   free(lineBuffers);
   upc_barrier;
   
   if (!MYTHREAD)
      printf("Overall time for %s is %.2f s\n", basename(argv[0]),
             ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
   
   return 0;
}
