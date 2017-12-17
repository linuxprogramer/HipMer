#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "../../common/optlist.h"
#include <upc.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <libgen.h>

#include "../../common/upc_compatibility.h"
#include "../../common/common.h"
#include "spannerUtils.h"

shared int64_t singletons = 0;
shared int64_t fiveRejects = 0;
shared int64_t threeRejects = 0;
shared int64_t uninformatives = 0;

int main(int argc, char **argv) {
    upc_tick_t start_time = upc_ticks_now();

   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "l:i:s:r:m:S:ARF:T:U:D:C:Z:B:");
   int reverseComplement = 0, innieRemoval = 0, minEndSeparation = 0, minFreqReported = 10, fivePrimeWiggleRoom = 5, threePrimeWiggleRoom = 5, truncate = 0, minMatch = 0, insertSize = 0, insertSigma = 0, readLength = 0, moreAlignmentsExist;
   char outputMerAligner[MAX_FILE_PATH];
   char *libname = NULL;
   FILE *laneFD1, *laneFD2, *srfFD, *outFD;
   int64_t i, j, k;
   align_info result1, result2;
   char *lineBuffers = (char*) calloc_chk(MAX_LINE_SIZE * 12, 1);
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
   char *srfLine = newReadIdLane2 + MAX_LINE_SIZE;
#ifdef MERGE
   char *aux = srfLine + MAX_LINE_SIZE;
   FILE *mergedFD = NULL;
#endif
   char *resRead1, *resRead2;
   int validAlign1, validAlign2;
   char *srfSuffixName = NULL;
   char srfFileName[MAX_FILE_PATH];
   UPC_TICK_T start, end;
   int shortPair = 600, trimOff;
   int subjectID;
   int count, endDistance;
   int64_t totalContigs = 0, totalScaffolds = 0, pairsProcessed;
   shared[1] contigScaffoldMap_t *contigScaffoldMap;
   contigScaffoldMap_t contigScaffoldMapEntry;
   shared[1] int64_t *scaffLengths;
   int splitRes, scaffoldId, contigId, cStrand, sStart, sEnd;
   int srfInput = 0;
   int subject1, projectedEnd1, projectedStart1, startStatus1, endStatus1, location1, simplePos, evenPos, oddPos;
   int qStart1, qStop1, qLength1, sStart1, sStop1, sLength1, strand1;
   int qStart2, qStop2, qLength2, sStart2, sStop2, sLength2, strand2;
   int subject2, projectedEnd2, projectedStart2, startStatus2, endStatus2, location2;
   char outputSpanner[MAX_FILE_PATH];
   int64_t singleton = 0, uninformative = 0, threeReject = 0, fiveReject = 0;
   
   double read_io_time = 0.0;
   UPC_TICK_T start_io, end_io;

   const char *base_dir = ".";
   
#ifdef MERGE
   {
       char myname[MAX_FILE_PATH];
       sprintf(myname, "NEWmerged_%d", MYTHREAD);
       get_rank_path(myname, MYTHREAD);
       mergedFD = fopen_chk(myname, "w");
   }
#endif
   
   
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'l':
            libname = thisOpt->argument;
            break;
         case 'S':
            /* FIXME: Standardize srf file output name conventions */
            srfInput = 1;
            srfSuffixName = thisOpt->argument;
            break;
         case 'U':
            truncate = atoi(thisOpt->argument);
            break;
         case 'T':
            threePrimeWiggleRoom = atoi(thisOpt->argument);
            break;
         case 'D':
            minEndSeparation = atoi(thisOpt->argument);
            break;
         case 'F':
            fivePrimeWiggleRoom = atoi(thisOpt->argument);
            break;
         case 'm':
            minMatch = atoi(thisOpt->argument);
            break;
         case 'i':
            insertSize = atoi(thisOpt->argument);
            break;
         case 's':
            insertSigma = atoi(thisOpt->argument);
            break;
         case 'r':
            readLength = atoi(thisOpt->argument);
            break;
         case 'A':
            innieRemoval = 1;
            break;
         case 'C':
            totalContigs = __atoi64(thisOpt->argument);
            break;
         case 'Z':
            totalScaffolds = __atoi64(thisOpt->argument);
            break;
         case 'R':
            reverseComplement = 1;
            break;
         case 'B':
             base_dir = thisOpt->argument;
             break;
         default:
            break;
      }
      free(thisOpt);
   }
   if (libname == NULL) 
       SDIE("You must specify a -l libname!\n");
   if (srfInput) {
       if (!srfSuffixName)
           SDIE("You must specify a -S srfSuffixName\n");
       sprintf(srfFileName, "%s/%s_%d", base_dir, srfSuffixName, MYTHREAD);
       get_rank_path(srfFileName, MYTHREAD);
       srfFD = fopen_chk(srfFileName, "r");
       if (!totalScaffolds || !totalContigs)
           SDIE("You must specify a -C totalContigs and a -Z totalScaffolds\n");
   }
   sprintf(outputMerAligner, "%s/%s-merAlignerOutput_%d_Read1", base_dir, libname, MYTHREAD);
   get_rank_path(outputMerAligner, MYTHREAD);
   laneFD1 = fopen_chk(outputMerAligner, "r");

   sprintf(outputMerAligner, "%s/%s-merAlignerOutput_%d_Read2", base_dir, libname, MYTHREAD);
   get_rank_path(outputMerAligner, MYTHREAD);
   laneFD2 = fopen_chk(outputMerAligner, "r");

   sprintf(outputSpanner, "%s/%s-spans_%d", base_dir, libname, MYTHREAD);
   get_rank_path(outputSpanner, MYTHREAD);
   outFD = fopen_chk(outputSpanner, "w+");
   
   endDistance = insertSize + 3 * insertSigma;
   char *fgets_result;
   upc_barrier;
   start = UPC_TICKS_NOW();
   
   /******************************/
   /*  Read scaffold report file */
   /******************************/
   if (srfInput) {
      /* Build scaffLengths data structure */
      if(MYTHREAD==0) printf("all allocating %ld + %ld\n", totalScaffolds*sizeof(int64_t), totalContigs*sizeof(contigScaffoldMap_t));
      scaffLengths = (shared[1] int64_t*) upc_all_alloc(totalScaffolds+1, sizeof(int64_t));
      if (scaffLengths == NULL) 
          DIE("Could not allocate %ld scaffolds!\n", totalScaffolds);
      for (i = MYTHREAD; i < totalScaffolds; i += THREADS) {
         scaffLengths[i] = 0;
      }
      
      /* Build contigScaffoldMap data structure */
      contigScaffoldMap = (shared[1] contigScaffoldMap_t*) upc_all_alloc(totalContigs+1,sizeof(contigScaffoldMap_t));
      if (contigScaffoldMap == NULL) 
          DIE("Could not allocate %ld totalContigs!\n", totalContigs); 
      for (i = MYTHREAD; i < totalContigs; i += THREADS) {
         contigScaffoldMap[i].scaffID = UNDEFINED;
      }
      
      upc_barrier;
      
      start_io = UPC_TICKS_NOW();
      fgets_result = fgets(srfLine, MAX_LINE_SIZE, srfFD);
      end_io = UPC_TICKS_NOW();
      read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
      
      while ( fgets_result  != NULL) {
         assert( fgets_result[ strlen(fgets_result) - 1 ] == '\n' );
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
         
         
         start_io = UPC_TICKS_NOW();

         fgets_result = fgets(srfLine, MAX_LINE_SIZE, srfFD);
         
         end_io = UPC_TICKS_NOW();
         read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
      }
      
      fclose(srfFD);
   }
   
   upc_barrier;
   
   /***************************/
   /* Read merAligner outputs */
   /***************************/
   start_io = UPC_TICKS_NOW();

   resRead1 = fgets(firstAlignmentBuffer1, MAX_LINE_SIZE, laneFD1);
   resRead2 = fgets(firstAlignmentBuffer2, MAX_LINE_SIZE, laneFD2);
   
   end_io = UPC_TICKS_NOW();
   read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
   moreAlignmentsExist = 1;

   while ( moreAlignmentsExist ) {
      
      if ((resRead1 == NULL) && (resRead2 == NULL)) break;
      
      /* Parse the first alignment of read (lane 1) */
      if (resRead1 != NULL) {
         assert( resRead1[ strlen(resRead1) - 1] == '\n');
         /* Split alignment and check for guard values */
#ifdef MERGE
         fprintf(mergedFD, "%s", firstAlignmentBuffer1);
#endif
         
         validAlign1 = splitAlignment(firstAlignmentBuffer1, readIdLane1, &qStart1, &qStop1, &qLength1, &subject1, &sStart1, &sStop1, &sLength1, &strand1);
         
         
         /* Assess alignment for completeness (do this before scaffold coordinate conversion!) */
         if (validAlign1 == SUCCESS) {
            validAlign1 = assessAlignment(strand1, qStart1, qStop1, qLength1, sStart1, sStop1, sLength1, fivePrimeWiggleRoom, threePrimeWiggleRoom, truncate, &projectedStart1, &projectedEnd1, &startStatus1, &endStatus1, &fiveReject, &threeReject);
         }
         
         /* Reorient alignment if requested */
         if (validAlign1 == SUCCESS) {
            if (reverseComplement) {
               strand1 = (strand1 == PLUS) ? MINUS : PLUS ;
               qStart1 = qLength1 - qStop1+1;
               qStop1 = qLength1 - qStart1+1;
            }
         }
         
         /* Convert to scaffold coordinate system if srfFile is specified */
         if (validAlign1 == SUCCESS) {
            if (srfInput) {
               validAlign1 = convertToScaffoldCoordinates(qStart1, qStop1, qLength1, &subject1, &strand1, &sStart1, &sStop1, &sLength1, &projectedStart1, &projectedEnd1, contigScaffoldMap, scaffLengths);
            }
         }
         
         /* Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle) */
         if (validAlign1 == SUCCESS) {
            validAlign1 = assessLocationAlignment(strand1, projectedStart1, projectedEnd1, sLength1, endDistance, &location1, &uninformative);
         }
         
         /* Store alignment in data structure result1 */
         if (validAlign1 == SUCCESS) {
            result1.subjects_matched = 0;
            simplePos = result1.subjects_matched;
            evenPos = 2 * simplePos;
            oddPos = evenPos + 1;
            result1.query_name = readIdLane1;
            result1.query_length = qLength1;
            result1.aligned_query_locs[evenPos] = qStart1;
            result1.aligned_query_locs[oddPos] = qStop1;
            result1.aligned_subject_ids[simplePos] = subject1;
            result1.aligned_subject_locs[evenPos] = sStart1;
            result1.aligned_subject_locs[oddPos] = sStop1;
            result1.aligned_strand[simplePos] = strand1;
            result1.startStatus[simplePos] = startStatus1;
            result1.endStatus[simplePos] = endStatus1;
            result1.location[simplePos] = location1;
            result1.aligned_subject_lengths[simplePos] = sLength1;
            result1.subjects_matched = 1;
         } else {
            result1.subjects_matched = 0;
         }
      }

      /* Read the the rest alignments of read (lane 1) */
      if (resRead1 != NULL) {
         
         start_io = UPC_TICKS_NOW();
         resRead1 = fgets(alignmentBuffer1, MAX_LINE_SIZE, laneFD1);
         end_io = UPC_TICKS_NOW();
         read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
         
         memcpy(copyAlignmentBuffer1, alignmentBuffer1, MAX_LINE_SIZE * sizeof(char));
         if (resRead1 != NULL) {
            assert( resRead1[ strlen(resRead1) - 1] == '\n');
           
#ifdef MERGE
            memcpy(aux, copyAlignmentBuffer1, MAX_LINE_SIZE * sizeof(char) );
#endif
            /* Split alignment and check for guard values */
            validAlign1 = splitAlignment(copyAlignmentBuffer1, newReadIdLane1, &qStart1, &qStop1, &qLength1, &subject1, &sStart1, &sStop1, &sLength1, &strand1);
            
            while ( strcmp(readIdLane1, newReadIdLane1) == 0 ) {
               
#ifdef MERGE
               fprintf(mergedFD, "%s", aux);
#endif
               
               /* Assess alignment for completeness (do this before scaffold coordinate conversion!) */
               if (validAlign1 == SUCCESS) {
                  validAlign1 = assessAlignment(strand1, qStart1, qStop1, qLength1, sStart1, sStop1, sLength1, fivePrimeWiggleRoom, threePrimeWiggleRoom, truncate, &projectedStart1, &projectedEnd1, &startStatus1, &endStatus1, &fiveReject, &threeReject);
               }
               
               /* Reorient alignment if requested */
               if (validAlign1 == SUCCESS) {
                  if (reverseComplement) {
                     strand1 = (strand1 == PLUS) ? MINUS : PLUS ;
                     qStart1 = qLength1 - qStop1+1;
                     qStop1 = qLength1 - qStart1+1;
                  }
               }
              
               /* Convert to scaffold coordinate system if srfFile is specified */
               if (validAlign1 == SUCCESS) {
                  if (srfInput) {
                     validAlign1 = convertToScaffoldCoordinates(qStart1, qStop1, qLength1, &subject1, &strand1, &sStart1, &sStop1, &sLength1, &projectedStart1, &projectedEnd1, contigScaffoldMap, scaffLengths);
                  }
               }
               
               /* Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle) */
               if (validAlign1 == SUCCESS) {
                  validAlign1 = assessLocationAlignment(strand1, projectedStart1, projectedEnd1, sLength1, endDistance, &location1, &uninformative);
               }
               
               /* Store alignment in data structure result1 */
               if (validAlign1 == SUCCESS) {
                  simplePos = result1.subjects_matched;
                  evenPos = 2 * simplePos;
                  oddPos = evenPos + 1;
                  result1.query_name = readIdLane1;
                  result1.query_length = qLength1;
                  result1.aligned_query_locs[evenPos] = qStart1;
                  result1.aligned_query_locs[oddPos] = qStop1;
                  result1.aligned_subject_ids[simplePos] = subject1;
                  result1.aligned_subject_locs[evenPos] = sStart1;
                  result1.aligned_subject_locs[oddPos] = sStop1;
                  result1.aligned_strand[simplePos] = strand1;
                  result1.aligned_subject_lengths[simplePos] = sLength1;
                  result1.startStatus[simplePos] = startStatus1;
                  result1.endStatus[simplePos] = endStatus1;
                  result1.location[simplePos] = location1;
                  result1.subjects_matched += 1;
               }
               
               start_io = UPC_TICKS_NOW();

               resRead1 = fgets(alignmentBuffer1, MAX_LINE_SIZE, laneFD1);
               
               end_io = UPC_TICKS_NOW();
               read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
               
               memcpy(copyAlignmentBuffer1, alignmentBuffer1, MAX_LINE_SIZE * sizeof(char));
               if (resRead1 == NULL) break;
               assert( resRead1[ strlen(resRead1) - 1] == '\n');
               validAlign1 = splitAlignment(copyAlignmentBuffer1, newReadIdLane1, &qStart1, &qStop1, &qLength1, &subject1, &sStart1, &sStop1, &sLength1, &strand1);
            }
         
            if (resRead1 != NULL) {
               memcpy(firstAlignmentBuffer1, alignmentBuffer1, MAX_LINE_SIZE * sizeof(char) );
               //strcpy(readIdLane1, newReadIdLane1);
            }
         }
      }
      
      /* Parse the first alignment of read (lane 2) */
      if (resRead2 != NULL) {
         assert( resRead2[ strlen(resRead2) - 1] == '\n');
         /* Split alignment and check for guard values */
#ifdef MERGE
         fprintf(mergedFD, "%s", firstAlignmentBuffer2);
#endif
         validAlign2 = splitAlignment(firstAlignmentBuffer2, readIdLane2, &qStart2, &qStop2, &qLength2, &subject2, &sStart2, &sStop2, &sLength2, &strand2);
         
         /* Assess alignment for completeness (do this before scaffold coordinate conversion!) */
         if (validAlign2 == SUCCESS) {
            validAlign2 = assessAlignment(strand2, qStart2, qStop2, qLength2, sStart2, sStop2, sLength2, fivePrimeWiggleRoom, threePrimeWiggleRoom, truncate, &projectedStart2, &projectedEnd2, &startStatus2, &endStatus2, &fiveReject, &threeReject);
         }
         
         /* Reorient alignment if requested */
         if (validAlign2 == SUCCESS) {
            if (reverseComplement) {
               strand2 = (strand2 == PLUS) ? MINUS : PLUS ;
               qStart2 = qLength2 - qStop2+1;
               qStop2 = qLength2 - qStart2+1;
            }
         }
         
         /* Convert to scaffold coordinate system if srfFile is specified */
         if (validAlign2 == SUCCESS) {
            if (srfInput) {
               validAlign2 = convertToScaffoldCoordinates(qStart2, qStop2, qLength2, &subject2, &strand2, &sStart2, &sStop2, &sLength2, &projectedStart2, &projectedEnd2, contigScaffoldMap, scaffLengths);
            }
         }
         
         /* Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle) */
         if (validAlign2 == SUCCESS) {
            validAlign2 = assessLocationAlignment(strand2, projectedStart2, projectedEnd2, sLength2, endDistance, &location2, &uninformative);
         }
         
         /* Store alignment in data structure result1 */
         if (validAlign2 == SUCCESS) {
            result2.subjects_matched = 0;
            simplePos = result2.subjects_matched;
            evenPos = 2 * simplePos;
            oddPos = evenPos + 1;
            result2.query_name = readIdLane2;
            result2.query_length = qLength2;
            result2.aligned_query_locs[evenPos] = qStart2;
            result2.aligned_query_locs[oddPos] = qStop2;
            result2.aligned_subject_ids[simplePos] = subject2;
            result2.aligned_subject_locs[evenPos] = sStart2;
            result2.aligned_subject_locs[oddPos] = sStop2;
            result2.aligned_strand[simplePos] = strand2;
            result2.startStatus[simplePos] = startStatus2;
            result2.endStatus[simplePos] = endStatus2;
            result2.location[simplePos] = location2;
            result2.aligned_subject_lengths[simplePos] = sLength2;
            result2.subjects_matched = 1;
         } else {
            result2.subjects_matched = 0;
         }
      }
      
      /* Read the the rest alignments of read (lane 2) */
      if (resRead2 != NULL) {
         
         start_io = UPC_TICKS_NOW();

         resRead2 = fgets(alignmentBuffer2, MAX_LINE_SIZE, laneFD2);
         
         end_io = UPC_TICKS_NOW();
         read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
         
         memcpy(copyAlignmentBuffer2, alignmentBuffer2, MAX_LINE_SIZE * sizeof(char));
         if (resRead2 != NULL) {
            assert( resRead2[ strlen(resRead2) - 1] == '\n');
            
#ifdef MERGE
            memcpy(aux, copyAlignmentBuffer2, MAX_LINE_SIZE * sizeof(char) );
#endif
            /* Split alignment and check for guard values */
            validAlign2 = splitAlignment(copyAlignmentBuffer2, newReadIdLane2, &qStart2, &qStop2, &qLength2, &subject2, &sStart2, &sStop2, &sLength2, &strand2);
            
            while ( strcmp(readIdLane2, newReadIdLane2) == 0 ) {
               
#ifdef MERGE
               fprintf(mergedFD, "%s", aux);
#endif
               
               /* Assess alignment for completeness (do this before scaffold coordinate conversion!) */
               if (validAlign2 == SUCCESS) {
                  validAlign2 = assessAlignment(strand2, qStart2, qStop2, qLength2, sStart2, sStop2, sLength2, fivePrimeWiggleRoom, threePrimeWiggleRoom, truncate, &projectedStart2, &projectedEnd2, &startStatus2, &endStatus2, &fiveReject, &threeReject);
               }
               
               /* Reorient alignment if requested */
               if (validAlign2 == SUCCESS) {
                  if (reverseComplement) {
                     strand2 = (strand2 == PLUS) ? MINUS : PLUS ;
                     qStart2 = qLength2 - qStop2+1;
                     qStop2 = qLength2 - qStart2+1;
                  }
               }
               
               /* Convert to scaffold coordinate system if srfFile is specified */
               if (validAlign2 == SUCCESS) {
                  if (srfInput) {
                     validAlign2 = convertToScaffoldCoordinates(qStart2, qStop2, qLength2, &subject2, &strand2, &sStart2, &sStop2, &sLength2, &projectedStart2, &projectedEnd2, contigScaffoldMap, scaffLengths);
                  }
               }
               
               /* Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle) */
               if (validAlign2 == SUCCESS) {
                  validAlign2 = assessLocationAlignment(strand2, projectedStart2, projectedEnd2, sLength2, endDistance, &location2, &uninformative);
               }
               
               /* Store alignment in data structure result1 */
               if (validAlign2 == SUCCESS) {
                  simplePos = result2.subjects_matched;
                  evenPos = 2 * simplePos;
                  oddPos = evenPos + 1;
                  result2.query_name = readIdLane2;
                  result2.query_length = qLength2;
                  result2.aligned_query_locs[evenPos] = qStart2;
                  result2.aligned_query_locs[oddPos] = qStop2;
                  result2.aligned_subject_ids[simplePos] = subject2;
                  result2.aligned_subject_locs[evenPos] = sStart2;
                  result2.aligned_subject_locs[oddPos] = sStop2;
                  result2.aligned_strand[simplePos] = strand2;
                  result2.startStatus[simplePos] = startStatus2;
                  result2.endStatus[simplePos] = endStatus2;
                  result2.location[simplePos] = location2;
                  result2.aligned_subject_lengths[simplePos] = sLength2;
                  result2.subjects_matched += 1;
               }
               
               start_io = UPC_TICKS_NOW();
               
               resRead2 = fgets(alignmentBuffer2, MAX_LINE_SIZE, laneFD2);
               
               end_io = UPC_TICKS_NOW();
               read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
               
               memcpy(copyAlignmentBuffer2, alignmentBuffer2, MAX_LINE_SIZE * sizeof(char));
               if (resRead2 == NULL) break;
               assert( resRead2[ strlen(resRead2) - 1] == '\n');

               validAlign2 = splitAlignment(copyAlignmentBuffer2, newReadIdLane2, &qStart2, &qStop2, &qLength2, &subject2, &sStart2, &sStop2, &sLength2, &strand2);
            }
            
            if (resRead2 != NULL) {
               memcpy(firstAlignmentBuffer2, alignmentBuffer2, MAX_LINE_SIZE * sizeof(char) );
               //strcpy(readIdLane2, newReadIdLane2);
            }
         }
      }
      
      /* Process pair */
      pairsProcessed += processPair(&result1, &result2, outFD, innieRemoval, minEndSeparation, srfInput, &uninformative, &singleton);
      strcpy(readIdLane1, newReadIdLane1);
      strcpy(readIdLane2, newReadIdLane2);
   
   }
   
   upc_fence;
   fclose(outFD);
   
#ifdef MERGE
   fclose(mergedFD);
#endif
   UPC_ATOMIC_FADD_I64(&singletons, singleton);
   UPC_ATOMIC_FADD_I64(&fiveRejects, fiveReject);
   UPC_ATOMIC_FADD_I64(&threeRejects, threeReject);
   UPC_ATOMIC_FADD_I64(&uninformatives, uninformative);

   upc_barrier;
   end = UPC_TICKS_NOW();
   if (MYTHREAD == 0) {
      printf("Unused alignments:\n");
      printf("UNINFORMATIVE: %ld\n", uninformatives);
      printf("5-TRUNCATED: %ld\n", fiveRejects);
      printf("3-TRUNCATED: %ld\n", threeRejects);
      printf("SINGLETON: %ld\n", singletons);
      printf("\nTime for spanner : %d seconds\n", ((int)UPC_TICKS_TO_SECS(end-start)));
      
      char countFileName[MAX_FILE_PATH];
      sprintf(countFileName,"%s-bmaMeta-spans", libname);
      get_rank_path(countFileName, -1);
      FILE *countFD = fopen_chk(countFileName, "w+" );
      fprintf(countFD, "%d\t%d\t%d\t%d\t%s-spans\n", insertSize, insertSigma, innieRemoval, readLength, libname);
      
      printf("\nRead I/O time is : %f seconds\n", read_io_time);
      printf("Write I/O time is : %f seconds\n", write_io_time);
      printf("Pure computation time is : %f seconds\n", (((int)UPC_TICKS_TO_SECS(end-start))/1000000.0)-(write_io_time + read_io_time));


      fclose(countFD);
   }
   
   upc_barrier;

	if (!MYTHREAD)
		printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
			   ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
   return 0;
}
