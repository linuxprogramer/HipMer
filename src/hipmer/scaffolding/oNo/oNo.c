#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
#include "../../common/optlist.h"
#include <upc.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <upc_nb.h>
#include <upc_tick.h>
#include <libgen.h>

#include "../../common/common.h"
#include "../../common/upc_compatibility.h"
#include "oNoHash.h"
#include "oNoUtils.h"

//#define DEBUG

#define AUX_STRINGS_SIZE 100
#define MAX_LINE_SIZE 1024
#define TIE_FAST_SIZE 50

shared int64_t total_scaff_entries = 0;
shared int64_t total_scaffs = 0;
shared int64_t posInCopy = 0;

//typedef FILE* fileptr;
//typedef shared[] fileptr* sharedFilePtr;

static int cmp_int(const void *a, const void *b)
{
    return (*(int*)a - *(int*)b);
}

static int calcN50(int *lens, int num)
{
    int tot_len = 0;
    for (int i = 0; i < num; i++) 
        tot_len += lens[i];
    qsort(lens, num, sizeof(int), cmp_int);
    // compute the N50
    long running_tot = 0;
    double idx = 0.5;
    for (int i = num - 1; i >= 0; i--) {
        running_tot += lens[i];
        if (running_tot >= idx * tot_len) {
            int nx = (int)(idx * 100);
            return lens[i];
        }
    }
    return 0;
}

int main(int argc, char **argv) {
   upc_tick_t start_time = upc_ticks_now();
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "m:l:p:s:c:L:n:C:o:B:N:");
   UPC_TICK_T start, endt;
   
   char *linkDataFilePrefix, linkDataFileName[MAX_FILE_PATH], *scaffoldReportFilePrefix, scaffoldReportFileName[MAX_FILE_PATH], *contigReportFilePrefix, contigReportFileName[MAX_FILE_PATH], linkMetaFile[MAX_FILE_PATH];
   FILE *scaffReportFD = NULL, *contigReportFD = NULL, *linkMetaFD = NULL, *curLinkFileFD = NULL;
   int merSize, pairThreshold = 2;
   scaf_entry_t **cachedEntries;
   
   char *lineBuffers = (char*) calloc_chk(MAX_LINE_SIZE*5, 1);
   char *line = lineBuffers;
   char *aux_line = line + MAX_LINE_SIZE;
   char *entries1 = aux_line + MAX_LINE_SIZE;
   char *entries2 = entries1 + MAX_LINE_SIZE;
   char *entries3 = entries2 + MAX_LINE_SIZE;

   int suspendable, maxInsertSize, LinkFiles = 0, length;
   libInfoType libInfo[MAX_NLIBS];
   int64_t nTotalObjects = 0, totalLinks = 0, totalContigs = 0;
   int srfFlag = 0, crfFlag = 0, scaffoldID;
   oNoObject localONOObject;
   scaf_report_t localScafReportObject;
   
   char *linkFileNames = calloc_chk(MAX_FILE_PATH * (MAX_NLIBS+1), 1);
   char *curLinkFile = linkFileNames + (MAX_FILE_PATH*MAX_NLIBS);
   
   shared[BS] int64_t *depthHist;
   shared[BS] int64_t *maxDepth;
   
   depthHist = (shared[BS] int64_t*) upc_all_alloc(MAX_HISTO_SIZE, sizeof(int64_t));
   maxDepth = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
   if (depthHist == NULL || maxDepth == NULL) 
       DIE("Could not allocate basic data structures!\n"); 
   
   int curNlibs;
   shared[BS] oNoObject *onoObjectLengths;
   shared[BS] sharedoNoObjectPtr *sharingScratchpad;
   sharedoNoObjectPtr *cachedSharingScratchpad = (sharedoNoObjectPtr *) malloc_chk(THREADS * sizeof(sharedoNoObjectPtr));
   
   shared[BS] sharedDoublePtr *sharingScratchpad2;
   sharedDoublePtr *cachedSharingScratchpad2 = (sharedDoublePtr *) malloc_chk(THREADS * sizeof(sharedDoublePtr));

   shared[BS] sharedScaffPtr *sharingScratchpad3;
   sharedScaffPtr *cachedSharingScratchpad3 = (sharedScaffPtr *) malloc_chk(THREADS * sizeof(sharedScaffPtr));

   
   //shared[BS] sharedScaffEntriesPtr *sharingScratchpad4;

   sharedoNoObjectPtr remoteLengthArray;
   sharedDoublePtr remoteLengthArray2;
   sharedScaffPtr remoteLengthArray3;
   sharedScaffEntriesPtr remoteLengthArray4;

   shared[BS] double *contigDepths, *scaffDepth;
   char *token, *aux, cur_type, cur_sign;
   int64_t contig, depth, depthInfoAvailable = 0, maxbases, peakDepth, report_ind, gap_id, previous_scaff_id, contig_length, cur_scaff_id, cur_id, cur_f1, cur_f2, cur_f3, cur_f4, cur_depth;
   int cur_f5;
   int64_t i, agg_depth, running_scaff_depth, running_real_bases, running_scaff_length, total_entries, nSpans, nSplints, mySpans, mySplints, totalTiesData;
   shared[BS] scaf_report_t *scaffReport;
   scaf_report_t *scaffReportLocal;
   scaf_entry_t current_scaff_report[MAX_SCAFF_REPORT_LENGTH];
   shared[] scaf_entry_t *new_ptr;
   oNolink_t *local_buffs;
   int64_t *local_index;
   oNoLink_heap_t link_heap;
   shared[BS] sharedTiePtr *tieSharedHeap;
   shared[BS] int64_t *tieSharedHeapPtr;
   char cur_typeS[30], cur_key[300];
   oNolink_t new_entry;
   int64_t hashval;
   shared[] tie_t *remoteTieHeap;
   shared[] tie_t *tmp;

   int64_t datum_heap_pos, link_heap_pos, entries_in_my_stack, localTiesPtr, splint_stack_ptr, span_stack_ptr, maxLikelihoodSplint, maxSplintCount, minUncertaintySpan, minSpanUncertainty, nUsedSpans, nUsedSplints, nAnomalousSpans, nAnomalousSplints, aux_ptr, remote_pos, posInEndTies, runningPtr;
   int cur_e1, cur_e2, cur_e1Ending, cur_e2Ending, gapEstimate, gapUncertainty, minimumGapSize, nSplintGap, splintGapEstimate, maxSplintError,splintAnomaly, n0, g0, u0, n, g, u, splintGapDiff, nSpanGap, spanGapEstimate, spanGapUncertainty, maxSpanZ, spanAnomaly, nGapLinks, valid_link, nGapL, gapE, gapU, cPos, cur_end;
   char cur_ending;
   double spanGapZ, gapZ;
   oNoDatum_t *datum_heap;
   oNolinkList_t *linkList_heap;
   oNolink_hash_table_t *local_hashtable;
   oNolink_t *my_stack;
   splints_stack_entry *splints;
   spans_stack_entry *spans;
   tie_t *endTies_local;
   oNolinkList_t *cur_bucket;
   oNoDatum_t *cur_data;
   dataTie_t *tiesDataHeap;
   endTies_t *endTies;
   tie_t *localTieHeap;
   oNoObject *localONOObjectLengths;
   dataTie_t *dataArray;
   splintDataTie_t *splintDataArray;

   oNoObject *sortedByLen;
   char *endMarks, endMark, *outputPrefix;
   double fdepth = 0.0, fcur_depth = 0.0, running_scaff_fdepth = 0.0;
   int64_t nInsertedSuspensions;
   nInsertedSuspensions = 0;

   const char *base_dir = ".";

   UPC_TICK_T start_timer, end_timer;
   UPC_TICK_T comm_start_timer, comm_end_timer;
   UPC_TICK_T serS, serE;
   UPC_TICK_T while1S, while1E;
   UPC_TICK_T while2S, while2E;
   UPC_TICK_T writeS, writeE;

   double while1Time = 0.0 , while2Time = 0.0 , writeTime = 0.0;

   int n50 = 0;
   int N_ACTIVE_THREADS = 3;

   /* Process the input arguments */
   /* TODO: Add support for blessed Link file */
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'l':
            strcpy(linkMetaFile, thisOpt->argument);
            linkMetaFD = fopen_rank_path(linkMetaFile, "r", -1);
            break;
         case 'm':
            merSize = atoi(thisOpt->argument);
            break;
         case 'L':
            totalLinks = __atoi64(thisOpt->argument) + 1;
            break;
         case 'n':
            nTotalObjects = __atoi64(thisOpt->argument) + 1;
            break;
         case 'C':
            totalContigs = __atoi64(thisOpt->argument) + 1;
            break;
         case 's':
            scaffoldReportFilePrefix = thisOpt->argument;
            srfFlag = 1;
            break;
         case 'c':
            contigReportFilePrefix = thisOpt->argument;
            crfFlag = 1;
            break;
         case 'p':
            pairThreshold = atoi(thisOpt->argument);
            break;
         case 'o':
            outputPrefix = thisOpt->argument;
            break;
         case 'B':
             base_dir = thisOpt->argument;
             break;
         case 'N':
            N_ACTIVE_THREADS = atoi(thisOpt->argument);
            break;
         default:
            break;
      }
      free(thisOpt);
   }

   if (srfFlag) {
       sprintf(scaffoldReportFileName, "%s/%s_%d", base_dir, scaffoldReportFilePrefix , MYTHREAD);
       scaffReportFD = fopen_rank_path(scaffoldReportFileName, "r", MYTHREAD);
   }

   if (crfFlag) {
       sprintf(contigReportFileName, "%s/%s_%d", base_dir, contigReportFilePrefix , MYTHREAD);
       contigReportFD = fopen_rank_path(contigReportFileName, "r", MYTHREAD);
   }

   if (crfFlag + srfFlag != 1) 
       DIE("Please one and only one of either: -c merDepthReport or -s scaffoldReportFile\n");
   
   upc_barrier;
   start = UPC_TICKS_NOW();
   
   /* TODO: Add support for blessed Link file */

   int nLinkFiles = 0;
   /* Read linkData file and store library information */
   while ( fgets(line, MAX_LINE_SIZE, linkMetaFD) != NULL ) {
      strcpy(aux_line, line);
      splitLinkMetaLine(aux_line, &libInfo[nLinkFiles]);
      serial_printf("Read lib %d: insert size %d, stddev %d, link file %s\n", nLinkFiles, 
                    libInfo[nLinkFiles].insertSize, libInfo[nLinkFiles].stdDev, 
                    libInfo[nLinkFiles].linkFileNamePrefix);
      nLinkFiles++;
   }
   fclose(linkMetaFD);
   
   /* WARNING: My implementation assumes that the libs are sorted based on the insert sizes */
   // This doesn't matter if there is only one input lib per oNo
   maxInsertSize = libInfo[nLinkFiles-1].insertSize;
   suspendable = maxInsertSize/2;

   start_timer = UPC_TICKS_NOW();
   int64_t my_nTotalObjects;

   if (crfFlag == 1) {
      nTotalObjects = nTotalObjects > totalContigs ? nTotalObjects : totalContigs;
      if(totalContigs == 0) totalContigs = nTotalObjects;
      /* All thread allocate onoObjectLengths shared array for contigs  */
#ifdef DEBUG
      if (MYTHREAD == 0) { printf("All_allocating %lld of %lld = %lld, plus locally %lld\n", (long long) nTotalObjects, (long long) sizeof(oNoObject)+sizeof(double), (long long) nTotalObjects * (sizeof(oNoObject)+sizeof(double)), (long long) nTotalObjects * sizeof(oNoObject)); }
#endif
      onoObjectLengths = (shared[BS] oNoObject*) upc_all_alloc(nTotalObjects, sizeof(oNoObject));
      sharingScratchpad = (shared[BS] sharedoNoObjectPtr*) upc_all_alloc(THREADS, sizeof(sharedoNoObjectPtr));
      sharingScratchpad2 = (shared[BS] sharedDoublePtr*) upc_all_alloc(THREADS, sizeof(sharedDoublePtr));
      contigDepths = (shared[BS] double*) upc_all_alloc(nTotalObjects, sizeof(double));
      if(onoObjectLengths == NULL || sharingScratchpad == NULL || sharingScratchpad2 == NULL || contigDepths == NULL) 
          DIE("Could not allocate %ld objects + misc\n", nTotalObjects); 

      if ( MYTHREAD < N_ACTIVE_THREADS) {
         my_nTotalObjects = (nTotalObjects+N_ACTIVE_THREADS-1)/N_ACTIVE_THREADS;
         sharingScratchpad[MYTHREAD] = (shared[] oNoObject*) upc_alloc(my_nTotalObjects * sizeof(oNoObject));
         //sharingScratchpad2[0] = (shared[] double*) upc_alloc(nTotalObjects * sizeof(double));
         if (sharingScratchpad[MYTHREAD] == NULL)
             DIE("Could not alloate %ld objects locally\n", nTotalObjects);
      }
      upc_barrier;
      remoteLengthArray = sharingScratchpad[MYTHREAD % N_ACTIVE_THREADS];
      //remoteLengthArray2 = sharingScratchpad2[0];
      
      /* Read the contig report file */
      int maxContigId = 0;
      while ( fgets(line, MAX_LINE_SIZE, contigReportFD) != NULL ) {
         token = strtok_r(line, "\t", &aux);
         assert(token != NULL);
         contig = atoi(token+6);       // Since the string should be ContigXXXXXXXX
         if (contig >= nTotalObjects) 
             DIE("too many contigs %ld in contigReport (expecting %ld): %s\n", contig, nTotalObjects, line);

         token = strtok_r(NULL, "\t", &aux);
         assert(token != NULL);
         length = atoi(token);         // Extract contig's length
         
         length = length + merSize -1; // Well, we store nKmers instead of the length in the merDepth file
         token = strtok_r(NULL, "\t", &aux);
         
         localONOObject.length = length;
         localONOObject.my_id = contig;

         onoObjectLengths[contig] = localONOObject;
         
         if (token != NULL) {
            fdepth = atof(token);       // Extract contig's depth (should be stored as int)
            char strmodal[20];
            sprintf(strmodal, "%0.f", fdepth);
            depth = (int) atof(strmodal);
            bupc_atomicI64_fetchadd_strict(&depthHist[depth], (int64_t)length);  // Increase appropriate entry in histogram
            contigDepths[contig] = fdepth;
            depthInfoAvailable = 1;
         }
         if (contig > maxContigId) maxContigId = contig;
      }
      fclose(contigReportFD);

      /* Thread 0 finds modal depth and propagates the info through shared array */
      upc_barrier;
      if (MYTHREAD == 0) {
         
         maxbases = 0;
         for (i=0; i<MAX_HISTO_SIZE; i++){
            agg_depth = depthHist[i];
            if (agg_depth > maxbases) {
               peakDepth = i;
               maxbases = agg_depth;
            }
         }
         for (i=0; i<THREADS; i++) maxDepth[i] = peakDepth;
//#ifdef DEBUG
         printf("Depth information available. Modal depth is %"PRId64".\n", peakDepth);
//#endif
      }
      upc_barrier;
      peakDepth = maxDepth[MYTHREAD];
      /* ensure all threads have the depth context */
      depthInfoAvailable = peakDepth > 0 ? 1 : 0; 
   }
   
   /* If an existing scaffold report exists, load scaffolding information */
   if (srfFlag == 1) {
      /* Allocate distributed data structure to store scafold reports */
#ifdef DEBUG
      if (MYTHREAD == 0) {
        printf("Input is SRF file\n");
        printf("All_allocating %ld, plus locally %ld\n", nTotalObjects*(sizeof(scaf_report_t)+sizeof(oNoObject)+sizeof(double))+totalContigs*sizeof(double), nTotalObjects * sizeof(oNoObject) + nTotalObjects * sizeof(scaf_report_t));
      }
#endif
      upc_barrier;
      scaffReport = (shared[BS] scaf_report_t*) upc_all_alloc(nTotalObjects, sizeof(scaf_report_t));
      onoObjectLengths = (shared[BS] oNoObject*) upc_all_alloc(nTotalObjects, sizeof(oNoObject));
      sharingScratchpad = (shared[BS] sharedoNoObjectPtr*) upc_all_alloc(THREADS, sizeof(sharedoNoObjectPtr));
      sharingScratchpad2 = (shared[BS] sharedDoublePtr*) upc_all_alloc(THREADS, sizeof(sharedDoublePtr));
      sharingScratchpad3 = (shared[BS] sharedScaffPtr*) upc_all_alloc(THREADS, sizeof(sharedScaffPtr));
      //sharingScratchpad4 = (shared[BS] sharedScaffEntriesPtr*) upc_all_alloc(THREADS, sizeof(sharedScaffEntriesPtr));
      contigDepths = (shared[BS] double*) upc_all_alloc(totalContigs, sizeof(double));
      scaffDepth = (shared[BS] double*) upc_all_alloc(nTotalObjects, sizeof(double));
      if(scaffReport == NULL || onoObjectLengths == NULL || sharingScratchpad == NULL || 
         sharingScratchpad2 == NULL || sharingScratchpad3 == NULL ||
         contigDepths == NULL || scaffDepth == NULL) 
          DIE("Could not allocate %ld and scratchpads\n", nTotalObjects*4);

      if ( MYTHREAD < N_ACTIVE_THREADS) {
         my_nTotalObjects = (nTotalObjects+N_ACTIVE_THREADS-1)/N_ACTIVE_THREADS;
         sharingScratchpad[MYTHREAD] = (shared[] oNoObject*) upc_alloc(my_nTotalObjects * sizeof(oNoObject));
         sharingScratchpad3[MYTHREAD] = (shared[] scaf_report_t*) upc_alloc(my_nTotalObjects * sizeof(scaf_report_t));
         if (sharingScratchpad[MYTHREAD] == NULL || sharingScratchpad3[MYTHREAD] == NULL)
             DIE("Could not allocate %ld (%ld bytes) on thread0 scratchpads\n", 
                 nTotalObjects, nTotalObjects * (sizeof(oNoObject) + sizeof(scaf_report_t)));
      }
      upc_barrier;
      remoteLengthArray = sharingScratchpad[MYTHREAD % N_ACTIVE_THREADS];
      //remoteLengthArray2 = sharingScratchpad2[0];
      remoteLengthArray3 = sharingScratchpad3[MYTHREAD % N_ACTIVE_THREADS];
      
      report_ind = 0;
      gap_id = 1;
      previous_scaff_id = UNDEFINED;
      contig_length = 0;
      running_scaff_fdepth = 0.0;
      running_real_bases = 0;
      running_scaff_length = 0;
      
#ifdef DEBUG
      printf("Thread %d: Done with allocation\n", MYTHREAD);
#endif
      /* Read the scaffold report file and store the info in the appropriate distributed table */
      while ( fgets(line, MAX_LINE_SIZE, scaffReportFD) != NULL ) {
         //printf("Thread %d: Read %s", MYTHREAD, line);
         token = strtok_r(line, "\t", &aux);
         assert(token != NULL);
         cur_scaff_id = atoi(token+8);    // Since the string should be ScaffoldXXXXXXX
         assert(cur_scaff_id >= 0);
         token = strtok_r(NULL, "\t", &aux);
         assert(token != NULL);
         //cur_type = ( (*token) == 'C' ) ? CONTIG : GAP ;  // Since column1 either CONTIGX or GAPX
         if ( (*token) == 'C' ) {
            cur_type = CONTIG;
         } else {
            cur_type = GAP;
         }
         
         token = strtok_r(NULL, "\t", &aux);
         assert(token != NULL);
         if (cur_type == CONTIG) {
            //cur_sign = ((*token) == '+') ? PLUS : MINUS ; // Since column2 is +ContigXXXX or -ContigXXXX
            if ((*token) == '+') {
               cur_sign = PLUS;
            } else {
               cur_sign = MINUS;
            }
            cur_id = atoi(token+7);                       // Since column2 is +ContigXXXX or -ContigXXXX
            token = strtok_r(NULL, "\t", &aux);
            assert(token != NULL);
            cur_f1 = atoi(token);
            token = strtok_r(NULL, "\t", &aux);
            assert(token != NULL);
            cur_f2 = atoi(token);
            token = strtok_r(NULL, "\t", &aux);
            // FIXME: Assumes that depth is always available
            if (token != NULL) {
               fcur_depth = atof(token);  // TODO: May want to read depth as double
               depthInfoAvailable = 1;
            } else {
               depthInfoAvailable = 0;
            }
            
         } else if (cur_type == GAP){
            //gap_id++;
            cur_f1 = atoi(token);
            token = strtok_r(NULL, "\t", &aux);
            assert(token != NULL);
            cur_f2 = atoi(token);
         }
         
         if ((cur_scaff_id == previous_scaff_id) || (previous_scaff_id == UNDEFINED)) {
            /* This is either the very first line of the scaffold report OR same line for PREVIOUS scaffold  */
            if (cur_type == CONTIG) {
               contig_length = cur_f2 - cur_f1 + 1;
               if ( cur_f2 > running_scaff_length) {
                  running_scaff_length = cur_f2;
               }

               if (depthInfoAvailable == 1) {
                  if (cur_id >= totalContigs) 
                      DIE("too many contigs (cur_id = %ld) in scaffold (expecting %ld): %s\n", cur_id, totalContigs, line);
                  contigDepths[cur_id] = fcur_depth;
                  running_real_bases += contig_length;
                  running_scaff_fdepth += contig_length * fcur_depth;
                  cur_depth = (int) trunc(fcur_depth);
                  char strmodal[20];
                  sprintf(strmodal, "%0.f", fcur_depth);
                  cur_depth = (int) atof(strmodal);
                  bupc_atomicI64_fetchadd_strict(&depthHist[cur_depth], (int64_t)contig_length);
               }
               
            }
            
            current_scaff_report[report_ind].type = cur_type;
            current_scaff_report[report_ind].f1 = cur_f1;
            current_scaff_report[report_ind].f2 = cur_f2;
            if (cur_type == CONTIG) {
               current_scaff_report[report_ind].entry_id = cur_id;
               //current_scaff_report[report_ind].intDepth = cur_depth;
               if (depthInfoAvailable == 1) {
                  current_scaff_report[report_ind].fDepth = fcur_depth;
               }
               current_scaff_report[report_ind].sign = cur_sign;
            } else if (cur_type == GAP) {
               current_scaff_report[report_ind].entry_id = gap_id;
               gap_id++;
            }
            report_ind++;
            previous_scaff_id = cur_scaff_id;
         } else {
            /* This is a line involving a NEW scaff report, store previous one in distributed data structure */
            new_ptr = (shared[] scaf_entry_t*) upc_alloc(report_ind * sizeof(scaf_entry_t));
            if (new_ptr == NULL) 
                DIE("Could not allocate %ld entries\n", report_ind);
            upc_memput(new_ptr, current_scaff_report, report_ind * sizeof(scaf_entry_t));
            UPC_ATOMIC_FADD_I64(&total_scaff_entries, report_ind );
            UPC_ATOMIC_FADD_I64(&total_scaffs, 1 );

            localScafReportObject.entries = new_ptr;
            localScafReportObject.nEntries = report_ind;
            localScafReportObject.realBases = running_real_bases;
            if (depthInfoAvailable == 1) {
               localScafReportObject.fdepth = running_scaff_fdepth;
            }
            localScafReportObject.length = running_scaff_length;
            if (previous_scaff_id >= nTotalObjects) 
                DIE("too many scaffolds  %ld (expecting %ld)\n", previous_scaff_id, nTotalObjects);
            scaffReport[previous_scaff_id] = localScafReportObject;
            
            localONOObject.my_id = previous_scaff_id;
            localONOObject.length = running_scaff_length;
            
            onoObjectLengths[previous_scaff_id] = localONOObject;
            
            //onoObjectLengths[previous_scaff_id].length = running_scaff_length;
            //contigDepths[previous_scaff_id] = running_scaff_fdepth/(1.0*running_real_bases);
            if (depthInfoAvailable == 1) {
               scaffDepth[previous_scaff_id] = running_scaff_fdepth/(1.0*running_real_bases);
            }

            //scaffReport[previous_scaff_id].entries = new_ptr;
            //scaffReport[previous_scaff_id].nEntries = report_ind;
            //scaffReport[previous_scaff_id].realBases = running_real_bases;
            //scaffReport[previous_scaff_id].depth = running_scaff_depth;
            //scaffReport[previous_scaff_id].length = running_scaff_length;

            /* Initialize again a current scaff report */
            running_real_bases = 0;
            running_scaff_depth = 0;
            running_scaff_fdepth = 0.0;
            running_scaff_length = 0;
            report_ind = 0;
            gap_id = 1;
            
            if (cur_type == CONTIG) {
               contig_length = cur_f2 - cur_f1 + 1;
               if ( cur_f2 > running_scaff_length) {
                  running_scaff_length = cur_f2;
               }
               
               if (depthInfoAvailable == 1) {
                  if (cur_id >= totalContigs) 
                      DIE("too many contigs in scaffoldreport (cur_id = %ld) (expecting %ld): %s\n", cur_id, totalContigs, line);
                  contigDepths[cur_id] = fcur_depth;
                  running_real_bases += contig_length;
                  running_scaff_fdepth += contig_length * fcur_depth;
                  cur_depth = (int) trunc(fcur_depth);
                  bupc_atomicI64_fetchadd_strict(&depthHist[cur_depth], (int64_t)contig_length);
               }
            }
            
            current_scaff_report[report_ind].type = cur_type;
            current_scaff_report[report_ind].f1 = cur_f1;
            current_scaff_report[report_ind].f2 = cur_f2;
            if (cur_type == CONTIG) {
               current_scaff_report[report_ind].entry_id = cur_id;
               //current_scaff_report[report_ind].intDepth = cur_depth;
               if (depthInfoAvailable == 1) {
                  current_scaff_report[report_ind].fDepth = fcur_depth;
               }
               
               current_scaff_report[report_ind].sign = cur_sign;
            } else if (cur_type == GAP) {
               current_scaff_report[report_ind].entry_id = gap_id;
               gap_id++;
            }
            report_ind++;
            previous_scaff_id = cur_scaff_id;
         }
      }
      fclose(scaffReportFD);

#ifdef DEBUG
      printf("Thread %d: Finished reading scaffoldReportFile. report_ind: %"PRId64"\n", MYTHREAD, report_ind);
#endif
      
      /* Check if the last scaffold report has been stored in the distributed scaffold report table */
      if (report_ind != 0) {
         new_ptr = (shared[] scaf_entry_t*) upc_alloc(report_ind * sizeof(scaf_entry_t));
         if (new_ptr == NULL) 
             DIE("Could not allocate %ld entries\n", report_ind);

         upc_memput(new_ptr, current_scaff_report, report_ind * sizeof(scaf_entry_t));
         UPC_ATOMIC_FADD_I64(&total_scaff_entries, report_ind );
         UPC_ATOMIC_FADD_I64(&total_scaffs, 1 );
         
         localScafReportObject.entries = new_ptr;
         localScafReportObject.nEntries = report_ind;
         localScafReportObject.realBases = running_real_bases;
         if (depthInfoAvailable == 1) {
            localScafReportObject.fdepth = running_scaff_fdepth;
         }
         localScafReportObject.length = running_scaff_length;
         if (previous_scaff_id >= nTotalObjects) 
             DIE("Too many scaffolds2 %ld (expecting %ld\n", previous_scaff_id, nTotalObjects);
         scaffReport[previous_scaff_id] = localScafReportObject;
         localONOObject.my_id = previous_scaff_id;
         localONOObject.length = running_scaff_length;
         onoObjectLengths[previous_scaff_id] = localONOObject;
         //onoObjectLengths[previous_scaff_id].length = running_scaff_length;
         if (depthInfoAvailable == 1) {
            scaffDepth[previous_scaff_id] = running_scaff_fdepth/(1.0*running_real_bases);
         }
         //scaffReport[previous_scaff_id].entries = new_ptr;
         //scaffReport[previous_scaff_id].nEntries = report_ind;
         //scaffReport[previous_scaff_id].realBases = running_real_bases;
         //running_real_bases = 0;
         //scaffReport[previous_scaff_id].depth = running_scaff_depth;
         //running_scaff_depth = 0;
         //scaffReport[previous_scaff_id].length = running_scaff_length;
         //onoObjectLengths[previous_scaff_id] = running_scaff_length;
         //running_scaff_length = 0;
      }
      
      /* Thread 0 finds modal depth and propagates the info through shared array */
      upc_barrier;
      
#ifdef DEBUG
      printf("Thread %d: Done with SRF file: %ld Scaffolds read\n", MYTHREAD, previous_scaff_id < 0 ? 0 : previous_scaff_id);
#endif
      /* Ask Jarrod about derived metric */
      /* Will compute this derived metric on the fly so no need to store it */
      /* Metric is contigDepths{} */
   }

   if (MYTHREAD == 0) {

      maxbases = 0;
      for (i=0; i<MAX_HISTO_SIZE; i++){
         agg_depth = depthHist[i];
         if (agg_depth > maxbases) {
            peakDepth = i;
            maxbases = agg_depth;
         }
      }
      for (i=0; i<THREADS; i++) maxDepth[i] = peakDepth;
#ifdef DEBUG
      printf("Thread %d: Depth information available. Modal depth is %ld.\n", MYTHREAD, peakDepth);
#endif
   }
   upc_barrier;
   peakDepth = maxDepth[MYTHREAD];
   /* ensure all threads have the depth context */
   depthInfoAvailable = peakDepth > 0 ? 1 : 0;

   
   /*  Build up contig linkage table using linkFiles */
   upc_barrier;
   
   end_timer = UPC_TICKS_NOW();
   double read_input_time = UPC_TICKS_TO_SECS(end_timer - start_timer);
   
   start_timer = UPC_TICKS_NOW();

   /* Allocate local buffers and local shared heaps required for the links hashtable */
   total_entries = totalLinks;
   /* Over-allocating in heaps to avoid OOM errors in case of large load imbalance */
   total_entries = (int64_t) total_entries * EXP_FACTOR;
   allocate_oNolocal_buffs(&local_buffs, &local_index);
   create_oNolink_heaps(total_entries, &link_heap);
   tieSharedHeap = (shared[BS] sharedTiePtr*) upc_all_alloc(THREADS, sizeof(sharedTiePtr));
   tieSharedHeapPtr = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
   upc_barrier;
   
   mySpans = nSpans = 0;
   mySplints = nSplints = 0;
   for ( i = 0; i < nLinkFiles; i++){
      sprintf(curLinkFile, "%s/%s_%d", base_dir, libInfo[i].linkFileNamePrefix, MYTHREAD);
      get_rank_path(curLinkFile, MYTHREAD);
#ifdef DEBUG
      printf("Reading linkFile: %s\n", curLinkFile);
#endif
      curLinkFileFD = fopen_chk( curLinkFile , "r");
      
      while ( fgets(line, MAX_LINE_SIZE, curLinkFileFD) != NULL ) {
          /* Parse the current link and store it */
         token = strtok_r(line, "\t", &aux);
         assert(token != NULL);
         strcpy(cur_typeS, token);
         token = strtok_r(NULL, "\t", &aux);
         assert(token != NULL);
         strcpy(cur_key, token);
         
         if ( cur_typeS[2] == 'A' ) {           // SP'A'N
            mySpans++;
            new_entry.link_type = SPAN;
         } else if ( cur_typeS[2] == 'L' ) {    // SP'L'INT
            mySplints++;
            new_entry.link_type = SPLINT;
         } else {
            printf("Thread %d: All links should be either spLint or spAn: %s\n", MYTHREAD, line);
            assert(0);
         }

         hashval = hashstr(THREADS, cur_key);
         new_entry.lib_id = i;
         
         token = strtok_r(NULL, "\t", &aux);
         assert(token != NULL);
         strcpy(entries1, token);
         token = strtok_r(NULL, "\t", &aux);
         assert(token != NULL);
         strcpy(entries2, token);
         
         if ( new_entry.link_type == SPAN ) {
            token = strtok_r(NULL, "\t", &aux);
            assert(token != NULL);
            strcpy(entries3, token);
         }
         
         if ( new_entry.link_type == SPAN ) {
            /* Link is of the form: SPAN	Scaffold10809495.3<=>Scaffold7516640.5	0|1	-42	0.87 */
            /* Parse first the key */
            token = strtok_r(cur_key, "<", &aux);
            assert(token != NULL);
            //new_entry.end1_s = ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) ? 3 : 5 ;
            if ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) {
               new_entry.end1_s = 3;
            } else {
               new_entry.end1_s = 5;
            }
            token[strlen(token)-2] = '\0';
            if ((*token) == 'S') {
               new_entry.end1_id = atoi(token + 8);
            } else {
               new_entry.end1_id = atoi(token + 6);
            }
            token = strtok_r(NULL, "<", &aux);
            assert(token != NULL);
            token += 2;    // to skip "=>"
            //new_entry.end2_s = ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) ? 3 : 5 ;
            if ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) {
               new_entry.end2_s = 3;
            } else {
               new_entry.end2_s = 5;
            }
            token[strlen(token)-2] = '\0';
            if ((*token) == 'S') {
               new_entry.end2_id = atoi(token + 8);
            } else {
               new_entry.end2_id = atoi(token + 6);
            }
            
            /* Parse then the info */
            token = strtok_r(entries1, "|", &aux);
            assert(token != NULL);
            new_entry.f1 = atoi(token);
            token = strtok_r(NULL, "|", &aux);
            assert(token != NULL);
            new_entry.f2 = atoi(token);
            new_entry.f3 = atoi(entries2);
            new_entry.f5 = atoi(entries3);
            
         } else if ( new_entry.link_type == SPLINT ) {
            /* Link is of the form: SPLINT	Contig21203848.5<=>Contig33389213.3	19|0|19	-46 */
            /* Parse first the key */
            token = strtok_r(cur_key, "<", &aux);
            assert(token != NULL);
            //new_entry.end1_s = ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) ? 3 : 5 ;
            if ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) {
               new_entry.end1_s = 3;
            } else {
               new_entry.end1_s = 5;
            }
            
            token[strlen(token)-2] = '\0';
            new_entry.end1_id = atoi(token + 6);
            token = strtok_r(NULL, "<", &aux);
            assert(token != NULL);
            token += 2;    // to skip "=>"
            //new_entry.end2_s = ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) ? 3 : 5 ;
            if ( strcmp( token+(strlen(token)-1) , "3" ) == 0 ) {
               new_entry.end2_s = 3;
            } else {
               new_entry.end2_s = 5;
            }
                
            token[strlen(token)-2] = '\0';
            new_entry.end2_id = atoi(token + 6);
            
            /* Parse then the info */
            token = strtok_r(entries1, "|", &aux);
            assert(token != NULL);
            new_entry.f1 = atoi(token);
            token = strtok_r(NULL, "|", &aux);
            assert(token != NULL);
            new_entry.f2 = atoi(token);
            token = strtok_r(NULL, "|", &aux);
            assert(token != NULL);
            new_entry.f3 = atoi(token);
            new_entry.f4 = atoi(entries2);
         }
         
         /* Store link to shared heaps */
         add_oNolink_to_shared_heaps(&new_entry, hashval, local_index, local_buffs, link_heap);
      }
      fclose(curLinkFileFD);
   
   }
   
   /* Adds remaining links to the shared heaps */
   add_rest_oNolinks_to_shared_heaps(local_index, local_buffs, link_heap);
   
   /* Wait until all heaps have been fixed */
   upc_barrier;
   
   end_timer = UPC_TICKS_NOW();
   double storing_links_time = UPC_TICKS_TO_SECS(end_timer - start_timer);
   
#ifdef DEBUG
   printf("Thread %d: Added links to shared heaps\n", MYTHREAD);
#endif
   //if (MYTHREAD == 0) {
   //   tieSharedHeap[0] = (shared[] tie_t*) upc_alloc(totalLinks * sizeof(tie_t));
   //   tieSharedHeapPtr[0] = 0;
   //}
   
   int64_t myProvision = (2 * totalLinks + THREADS -1)/THREADS + 1;
   tieSharedHeap[MYTHREAD] = (shared[] tie_t*) upc_alloc(myProvision * sizeof(tie_t));
   tieSharedHeapPtr[MYTHREAD] = 0;
   
   //remoteTieHeap = tieSharedHeap[0];
   upc_barrier;
   
   start_timer = UPC_TICKS_NOW();
   
   /* Create local hashtable by iterating on local heap and storing link_data to local_hashtable */
   datum_heap_pos = 0;
   link_heap_pos = 0;
   entries_in_my_stack = link_heap.heap_indices[MYTHREAD];
   create_oNolocal_heaps(entries_in_my_stack, &datum_heap, &linkList_heap);
   local_hashtable = create_oNolink_hash_table(entries_in_my_stack+1);
   my_stack = (oNolink_t*) link_heap.link_ptr[MYTHREAD];
   
#ifdef DEBUG
   printf("Thread %d: Entries in my stack: %ld\n", MYTHREAD, entries_in_my_stack);
   char myTies[MAX_FILE_PATH];
   sprintf(myTies, "myTies-p%d-%d.txt", pairThreshold, MYTHREAD);
   FILE *tiesFD = fopen_rank_path(myTies, "w", MYTHREAD);
#endif
   
   for (i=0; i < entries_in_my_stack; i++) {
      add_oNolink_data(local_hashtable, &my_stack[i], &datum_heap_pos, datum_heap, &link_heap_pos, linkList_heap);
   }
   
   
   end_timer = UPC_TICKS_NOW();
   double local_linktable_time = UPC_TICKS_TO_SECS(end_timer - start_timer);
   
#ifdef DEBUG
   printf("Thread %d: My splints: %ld, My spans: %ld\n", MYTHREAD, mySplints, mySpans);
   int64_t ntL = 0;
   double chS = 0.0;
#endif
   /* Allocate local stacks for splints and spans */
   nSpans = nSplints = entries_in_my_stack+1;
   splints = (splints_stack_entry*) malloc_chk(nSplints * sizeof(splints_stack_entry));
   spans = (spans_stack_entry*) malloc_chk(nSpans * sizeof(spans_stack_entry));
   
   /* Allocate a single local buffer for endTies(). All of them will go to processor 0 */
   endTies_local = (tie_t*) malloc_chk( MAX_TIES_LOCAL * sizeof(tie_t) );
   localTiesPtr = 0;
   
   int64_t totTies = 0;
   
   /**************************************************************************************************/
   /* Now iterate over the hash table and build ties between contigs consolidating splint/span links */
   /**************************************************************************************************/
   start_timer = UPC_TICKS_NOW();
   double build_ties_communication_time = 0.0;
   char thereIsASplint = 0;
   
   for (i=0; i < local_hashtable->size; i++) {
      cur_bucket = local_hashtable->table[i];
      while (cur_bucket != NULL) {
         
         //totTies++;

         
         /* Initialize splints and spans local stacks */
         splint_stack_ptr = 0;
         span_stack_ptr = 0;
         
         maxLikelihoodSplint = UNDEFINED;
         maxSplintCount = 0;
         minUncertaintySpan = UNDEFINED;
         minSpanUncertainty = 100000;
         
         /* TODO: Add support for Blessed link files */
         
         /* Process the data in the the current link */
         cur_e1 = cur_bucket->end1_id;
         cur_e2 = cur_bucket->end2_id;
         cur_e1Ending = cur_bucket->end1_s;
         cur_e2Ending = cur_bucket->end2_s;
         //cur_type = cur_bucket->link_type;
         cur_data = cur_bucket->data;
         thereIsASplint = 0;
         while (cur_data != NULL) {
            cur_type = cur_data->link_type;
            /* Different process of datum depending on its type (SPLINT/SPAN) */
            if ( cur_type == SPLINT ) {
               nUsedSplints = cur_data->f1;
               nAnomalousSplints = cur_data->f2;
               //nAllSplints = cur_data->f3;
               gapEstimate = cur_data->f4;
               //ntL++;
               thereIsASplint = 1;
               if ( nUsedSplints >= pairThreshold ) {
                  assert(splint_stack_ptr < nSplints);
                  assert(splint_stack_ptr >= 0);
                  splints[splint_stack_ptr].nUsedSplints = nUsedSplints;
                  splints[splint_stack_ptr].gapEstimate = gapEstimate;
                  splints[splint_stack_ptr].linkFileID = cur_data->lib_id;
                  if ( maxLikelihoodSplint == UNDEFINED ) {
                     maxLikelihoodSplint = splint_stack_ptr;
                     maxSplintCount = nUsedSplints;
                  }
                  if ( nUsedSplints > maxSplintCount ) {
                     maxLikelihoodSplint = splint_stack_ptr;
                     maxSplintCount = nUsedSplints;
                  }
                  splint_stack_ptr++;
               }
            } else if ( cur_type == SPAN ) {
               nAnomalousSpans = cur_data->f1;
               nUsedSpans = cur_data->f2;
               gapEstimate = cur_data->f3;
               gapUncertainty = cur_data->f5;
               //fprintf(tiesFD,"NL\n");
               //fprintf(tiesFD, "%ld|%ld,%d,%d\n",nAnomalousSpans ,nUsedSpans, gapEstimate, gapUncertainty );
               //ntL++;

               if ( nUsedSpans >= pairThreshold ) {
                  assert(span_stack_ptr < nSpans);
                  assert(span_stack_ptr >= 0);
                  spans[span_stack_ptr].nUsedSpans = nUsedSpans;
                  spans[span_stack_ptr].gapEstimate = gapEstimate;
                  spans[span_stack_ptr].gapUncertainty = gapUncertainty;
                  spans[span_stack_ptr].linkFileID = cur_data->lib_id;
                  if ( minUncertaintySpan == UNDEFINED ) {
                     minUncertaintySpan = span_stack_ptr;
                     minSpanUncertainty = gapUncertainty;
                  }
                  if ( gapUncertainty < minSpanUncertainty ) {
                     minUncertaintySpan = span_stack_ptr;
                     minSpanUncertainty = gapUncertainty;
                  }
                  span_stack_ptr++;
               }
            }
            cur_data = cur_data->next;
         }
         /* END: Process the data in the the current link */
         assert(nSplints >= splint_stack_ptr);
         assert(nSpans >= span_stack_ptr);
         
         minimumGapSize = -(merSize-2);
         nSplintGap = UNDEFINED;
         splintGapEstimate = UNDEFINED;
         maxSplintError = 2;
         splintAnomaly = 0;
         if ( maxLikelihoodSplint != UNDEFINED ) {
            n0 = splints[maxLikelihoodSplint].nUsedSplints;
            g0 = splints[maxLikelihoodSplint].gapEstimate;
            //splintGapEstimate = (g0 < minimumGapSize) ? minimumGapSize : g0 ;
            if ( g0 < minimumGapSize ) {
               splintGapEstimate = minimumGapSize;
            } else {
               splintGapEstimate = g0;
            }
            nSplintGap = n0;
            for (aux_ptr = 0; aux_ptr < splint_stack_ptr; aux_ptr++) {
               n = splints[aux_ptr].nUsedSplints;
               g = splints[aux_ptr].gapEstimate;
               splintGapDiff = abs(g-splintGapEstimate);
               if (splintGapDiff > maxSplintError) {
                  splintAnomaly = 1;
               }
            }
         }
         
         nSpanGap = UNDEFINED;
         spanGapEstimate = UNDEFINED;
         spanGapUncertainty = UNDEFINED;
         maxSpanZ = 3;

         spanAnomaly = 0;
         if ( minUncertaintySpan != UNDEFINED ) {
          //  totTies++;
            n0 = spans[minUncertaintySpan].nUsedSpans;
            g0 = spans[minUncertaintySpan].gapEstimate;
            u0 = spans[minUncertaintySpan].gapUncertainty;
            //fprintf(tiesFD,"minUncertaintySpan: %d, %d, %d \n", n0, g0, u0);

            //spanGapEstimate = (g0 < minimumGapSize) ? minimumGapSize : g0 ;
            if ( g0 < minimumGapSize ) {
               spanGapEstimate = minimumGapSize;
            } else {
               spanGapEstimate = g0;
            }
            //spanGapUncertainty = (u0 < 1) ? 1 : u0 ;
            if ( u0 < 1 ) {
               spanGapUncertainty = 1;
            } else {
               spanGapUncertainty = u0;
            }
            nSpanGap = n0;
            for (aux_ptr = 0; aux_ptr < span_stack_ptr; aux_ptr++) {
               n = spans[aux_ptr].nUsedSpans;
               g = spans[aux_ptr].gapEstimate;
               u = spans[aux_ptr].gapUncertainty;
               spanGapZ = abs(g - spanGapEstimate)/(1.0*(u+1));
               //chS += spanGapZ;
               if (spanGapZ > (maxSpanZ *1.0)) {
                  spanAnomaly = 1;
                  //totTies++;
               }
            }
         }
         
         gapEstimate = UNDEFINED;
         gapUncertainty = UNDEFINED;
         nGapLinks = UNDEFINED;
         valid_link = 1;

         /* if ( spanGapEstimate == UNDEFINED) {
            fprintf(tiesFD,"spanGapEstimate: \n");
         } else {
            fprintf(tiesFD,"spanGapEstimate: %d\n", spanGapEstimate);
         }*/

         
         if ( splintGapEstimate != UNDEFINED ) {
            gapEstimate = splintGapEstimate;
            gapUncertainty = 1;
            nGapLinks = nSplintGap;
            //totTies++;
         } else if ( spanGapEstimate != UNDEFINED) {
            gapEstimate = spanGapEstimate;
            gapUncertainty = spanGapUncertainty;
            nGapLinks = nSpanGap;
            //totTies++;
         } else {
            //totTies++;
            valid_link = 0;
         }
         
         if (valid_link == 1) {
            if ( (spanGapEstimate != UNDEFINED) && (splintGapEstimate != UNDEFINED) ) {
               gapZ = (1.0*(splintGapEstimate-spanGapEstimate))/(1.0*spanGapUncertainty);
               if ( fabs(gapZ) > maxSpanZ ) {
                  /* TODO: add warning message */
               }
            }
            
            endTies_local[localTiesPtr].end1_id = cur_e1;
            endTies_local[localTiesPtr].end2_id = cur_e2;
            endTies_local[localTiesPtr].ending1 = cur_e1Ending;
            endTies_local[localTiesPtr].ending2 = cur_e2Ending;
            endTies_local[localTiesPtr].nGapLinks = nGapLinks;
            endTies_local[localTiesPtr].gapEstimate = gapEstimate;
            endTies_local[localTiesPtr].gapUncertainty = gapUncertainty;
            endTies_local[localTiesPtr].splintFlag = thereIsASplint;
            localTiesPtr++;
            
#ifdef DEBUG
            //totTies++;
            //fprintf(tiesFD, "Tie: Contig%d.%d <=> Contig%d.%d | %d | %d | %d \n", cur_e1, (int) cur_e1Ending, cur_e2, (int) cur_e2Ending, nGapLinks, gapEstimate, gapUncertainty);
#endif
            
            if ( localTiesPtr == MAX_TIES_LOCAL ) {
               comm_start_timer = UPC_TICKS_NOW();
               
               //remote_pos = bupc_atomicI64_fetchadd_strict(&tieSharedHeapPtr[0], MAX_TIES_LOCAL);
               remote_pos = tieSharedHeapPtr[MYTHREAD];
               
               if (remote_pos + MAX_TIES_LOCAL >= myProvision) {
                  tmp = (shared[] tie_t*) upc_alloc(2 * (remote_pos + MAX_TIES_LOCAL) * sizeof(tie_t));
                  memcpy( (tie_t*) tmp, (tie_t*) tieSharedHeap[MYTHREAD], remote_pos * sizeof(tie_t));
                  upc_free(tieSharedHeap[MYTHREAD]);
                  tieSharedHeap[MYTHREAD] = tmp;
                  myProvision = 2 * (remote_pos + MAX_TIES_LOCAL);
               }
               
               memcpy( (tie_t*) (tieSharedHeap[MYTHREAD] + remote_pos), endTies_local, MAX_TIES_LOCAL * sizeof(tie_t) );
               tieSharedHeapPtr[MYTHREAD] += MAX_TIES_LOCAL;
               //upc_memput( remoteTieHeap + remote_pos, endTies_local, MAX_TIES_LOCAL * sizeof(tie_t) );
               
               localTiesPtr = 0;
               
               comm_end_timer = UPC_TICKS_NOW();
               build_ties_communication_time += UPC_TICKS_TO_SECS(comm_end_timer - comm_start_timer);
            }
            
            
            
         }
         cur_bucket = cur_bucket->next;
      }
   }
   
   /* Store remaining ties to remote heap */
   if ( localTiesPtr != 0 ) {
      comm_start_timer = UPC_TICKS_NOW();

      //remote_pos = bupc_atomicI64_fetchadd_strict(&tieSharedHeapPtr[0], localTiesPtr);
      //upc_memput( remoteTieHeap + remote_pos, endTies_local, localTiesPtr * sizeof(tie_t) );
      remote_pos = tieSharedHeapPtr[MYTHREAD];
      
      if (remote_pos + localTiesPtr >= myProvision) {
         tmp = (shared[] tie_t*) upc_alloc(2 * (remote_pos + localTiesPtr) * sizeof(tie_t));
         memcpy( (tie_t*) tmp, (tie_t*) tieSharedHeap[MYTHREAD], remote_pos * sizeof(tie_t));
         upc_free(tieSharedHeap[MYTHREAD]);
         tieSharedHeap[MYTHREAD] = tmp;
         myProvision = 2 * (remote_pos + localTiesPtr);
      }
      
      memcpy( (tie_t*) (tieSharedHeap[MYTHREAD] + remote_pos), endTies_local, localTiesPtr * sizeof(tie_t) );
      tieSharedHeapPtr[MYTHREAD] += localTiesPtr;
      
      comm_end_timer = UPC_TICKS_NOW();
      build_ties_communication_time += UPC_TICKS_TO_SECS(comm_end_timer - comm_start_timer);
   }
   
   if (splints != NULL) free(splints);
   if (spans != NULL) free(spans);
   
   end_timer = UPC_TICKS_NOW();
   double build_ties_time = UPC_TICKS_TO_SECS(end_timer - start_timer);
   
#ifdef DEBUG
   //fclose(tiesFD);
   printf("Thread %d: Checksum is %f. %ld total Links, %ld ties, %ld totalContigs %ld totalScaffolds\n", MYTHREAD, chS, ntL, tieSharedHeapPtr[0], totalContigs, nTotalObjects);
#endif
   
   start_timer = UPC_TICKS_NOW();
#ifdef DEBUG
   printf("Thread %d: past start_timer: depthInfoAvailable: %ld\n", MYTHREAD, depthInfoAvailable);
#endif
   
   upc_barrier;

#ifdef DEBUG
   printf("Thread %d: past variable: depthInfoAvailable: %ld\n", MYTHREAD, depthInfoAvailable);
#endif
   
   int64_t myTotalObjects;
   if (depthInfoAvailable == 1) {
      if (MYTHREAD < N_ACTIVE_THREADS) {
          //sharingScratchpad2[0] = (shared[] double*) upc_alloc(nTotalObjects * sizeof(double));
         myTotalObjects = (totalContigs + N_ACTIVE_THREADS -1) / N_ACTIVE_THREADS;
         sharingScratchpad2[MYTHREAD] = (shared[] double*) upc_alloc(myTotalObjects * sizeof(double));
         if (sharingScratchpad2[MYTHREAD] == NULL)
             DIE("Could not upc_alloc sharingScratchpad2 %ld doubles\n", nTotalObjects);
      }
   }
  
   // barrier must be on all threads, regardless if they set depthInfoAvailable...
   upc_barrier;

   if (depthInfoAvailable == 1) {
         if (sharingScratchpad2[MYTHREAD % N_ACTIVE_THREADS] == NULL)
             DIE("Thread %d: Could not use sharingScratchpad2 %ld doubles\n", MYTHREAD, nTotalObjects);
         //else 
         //    printf("Success!\n"); 
      
      remoteLengthArray2 = sharingScratchpad2[MYTHREAD % N_ACTIVE_THREADS];
   }

#ifdef DEBUG
   printf("Thread %d: Done allocating scratchpad\n", MYTHREAD);
#endif


   /* All threads store their oNoObjectlengths to thread 0 */
   //printf("nTotalObjects %ld\n", nTotalObjects);
   for (i=MYTHREAD; i<nTotalObjects; i+=THREADS) {

      //remoteLengthArray[i] = onoObjectLengths[i];
      remoteLengthArray[i/N_ACTIVE_THREADS] = onoObjectLengths[i];
      
      
      if (depthInfoAvailable == 1) {
         if (srfFlag) {
            //remoteLengthArray3[i] = scaffReport[i];
         } else {
            remoteLengthArray2[i/N_ACTIVE_THREADS] = contigDepths[i];
            //printf("Contig%ld %f\n", i, remoteLengthArray2[i]);
         }
      }
      //printf("Contig%ld length %d\n", i ,onoObjectLengths[i].length);
   }

   if (depthInfoAvailable == 1) {
      if (srfFlag == 1) {
          //printf("nTotalObjects %ld, contigs %ld\n", nTotalObjects, totalContigs);
          //for (i=MYTHREAD; i<nTotalObjects; i+=THREADS) {
          for (i=MYTHREAD; i<totalContigs; i+=THREADS) {
            remoteLengthArray2[i/N_ACTIVE_THREADS] = contigDepths[i];
            //printf("Contig%ld %f\n", i, remoteLengthArray2[i]);
         }
      }
   }

   upc_barrier;

#ifdef DEBUG
   printf("Thread %d: Done calculating depths\n", MYTHREAD);
#endif

   if (total_scaffs >= nTotalObjects) 
       SDIE("The total number of objects %ld is less than the number of scaffolds %ld\n"
           "Check the parameters: is -n set correctly?\n", nTotalObjects, total_scaffs);
   int64_t localBound = total_scaffs;
   int curEntries = 0;
   int64_t curPos = 0;
   
   if ( srfFlag == 1 ) {
      //if ( MYTHREAD == 0 ) {
      //   printf("Local allocating %ld of %lu\n", total_scaff_entries, (unsigned long)sizeof(scaf_entry_t));
      //   sharingScratchpad4[0] = (shared[] scaf_entry_t*) upc_alloc(total_scaff_entries * sizeof(scaf_entry_t));
      //}
      //upc_barrier;
      
      //remoteLengthArray4 = sharingScratchpad4[0];
      
      for ( i=MYTHREAD; i < localBound; i += THREADS ) {
         
         //curEntries = scaffReport[i].nEntries;
         //curPos = UPC_ATOMIC_FADD_I64(&posInCopy, curEntries );
         
         //if (curPos >= total_scaff_entries)
         //    DIE("too many scaff entries %ld (expecting %ld)\n", curPos, total_scaff_entries);
         //upc_memcpy(&remoteLengthArray4[curPos], scaffReport[i].entries, curEntries * sizeof(scaf_entry_t) );
         
         //scaffReport[i].entries = &remoteLengthArray4[curPos];
         remoteLengthArray3[i/N_ACTIVE_THREADS] = scaffReport[i];
         
      }
   
#ifdef DEBUG
   printf("Thread %d: Calculated the srf entries\n", MYTHREAD);
#endif
      upc_barrier;
   }
   
   scaf_report_t *tmpBuf;
   int64_t v, curEnt, cur;
   
   if ( (srfFlag == 1) && (MYTHREAD == 0) ) {
      scaffReportLocal = (scaf_report_t*) malloc_chk(localBound * sizeof(scaf_report_t));
      
      //for (i = 0; i< THREADS; i++) {
      //   cachedSharingScratchpad3[i] = sharingScratchpad3[i];
      //}
      
      curEnt = (nTotalObjects+N_ACTIVE_THREADS-1)/N_ACTIVE_THREADS ;
      tmpBuf = (scaf_report_t*) malloc_chk(curEnt * sizeof(scaf_report_t));
      
      for (v = 0; v < N_ACTIVE_THREADS; v++) {
         upc_memget(tmpBuf, sharingScratchpad3[v], curEnt * sizeof(scaf_report_t));
         cur = 0;
         for ( i = v; i < localBound; i+= N_ACTIVE_THREADS) {
            scaffReportLocal[i] = tmpBuf[cur];
            cur++;
         }
      
      }
      
      free(tmpBuf);
      
      //for (i = 0; i <localBound; i++) {
      //   scaffReportLocal[i] = *(cachedSharingScratchpad3[i%N_ACTIVE_THREADS] + i/N_ACTIVE_THREADS);
      //}
      
   }
   
   upc_barrier;

#ifdef DEBUG
   printf("Thread %d: Done with scf\n", MYTHREAD);
#endif
   
   end_timer = UPC_TICKS_NOW();
   double propagating_lengths_time = UPC_TICKS_TO_SECS(end_timer - start_timer);
   
   start_timer = UPC_TICKS_NOW();
   
   double sorting_time = 0.0;
   UPC_TICK_T start_sort, end_sort;
   double storTiesTime = 0.0, markendTime = 0.0, suspensionsTime = 0.0, bestTiesTime = 0.0, lockEndsTime = 0.0, traverseLocksTime = 0.0, fileOpeningTime = 0.0;
   
   double *localContigDepth = NULL;
   int64_t curEnt3, pos3;
   double *tmpBuf3;
   
   if (MYTHREAD == 0) {
       //assert(upc_threadof(remoteLengthArray2) == MYTHREAD);
      localContigDepth = (double*) malloc_chk(totalContigs * sizeof(double));
      
      //for (i = 0; i < THREADS; i++) {
      //   cachedSharingScratchpad2[i] = sharingScratchpad2[i];
      //}
      
      int64_t heapNumber;
      
      if (depthInfoAvailable == 1) {
         
         curEnt3 = (totalContigs + N_ACTIVE_THREADS -1) / N_ACTIVE_THREADS;
         tmpBuf3 = (double*) malloc_chk(curEnt3 * sizeof(double));
         
         for (heapNumber = 0; heapNumber < N_ACTIVE_THREADS; heapNumber++) {
            upc_memget(tmpBuf3, sharingScratchpad2[heapNumber], curEnt3 * sizeof(double));
            
            pos3 = 0;
            if (crfFlag) {
               for (i=heapNumber; i < nTotalObjects; i+= N_ACTIVE_THREADS ) {
                  localContigDepth[i] = tmpBuf3[pos3];
                  pos3++;
               }
            }
            
            if (srfFlag) {
               for (i=heapNumber; i < totalContigs; i+= N_ACTIVE_THREADS ) {
                  localContigDepth[i] = tmpBuf3[pos3];
                  pos3++;
               }
            }
            
         }
         
         free(tmpBuf3);
         
         //if (crfFlag) {
         //   for (i = 0; i < nTotalObjects; i++) {
         //      localContigDepth[i] = *(cachedSharingScratchpad2[i%N_ACTIVE_THREADS] + i/N_ACTIVE_THREADS);
         //   }
         //}
         //if (srfFlag) {
         //   for (i = 0; i < totalContigs; i++) {
         //      localContigDepth[i] = *(cachedSharingScratchpad2[i%N_ACTIVE_THREADS] + i/N_ACTIVE_THREADS);
         //   }
         //}
      }
   }

   /* Traverse end locks to build scaffolds */

#ifdef DEBUG
   printf("Thread %d: Allocating sharedScaffolds\n", MYTHREAD);
#endif
   
   /* Allocate shared buffers for storing the scaffold report files */
   shared[BS] ScaffoldReportFiles *sharedScaffolds = (shared[BS] ScaffoldReportFiles*) upc_all_alloc(THREADS, sizeof(ScaffoldReportFiles));
   
   //shared[BS] sharedCharPtr *sharedScaffoldsPointers = (shared[BS] sharedCharPtr*) upc_all_alloc(THREADS, sizeof(sharedCharPtr));
   //shared[BS] int64_t *sharedScaffoldsPositions = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
   //if (sharedScaffoldsPointers == NULL || sharedScaffoldsPositions == NULL) {
   if (sharedScaffolds == NULL) 
       DIE("Could not allocate sharedScaffoldsPointers and Positions\n");
   
#ifdef DEBUG
   printf("Thread %d: Allocated sharedScaffolds\n", MYTHREAD);
#endif
   
   /* Over-estimation of output size for each thread */
   int64_t lineSize = 400;
   int64_t expansion = 3;
   int64_t outputSize = (nTotalObjects * lineSize * expansion + THREADS - 1) / THREADS;
   sharedScaffolds[MYTHREAD].ptr = (sharedCharPtr) upc_alloc(outputSize * sizeof(char));
   sharedScaffolds[MYTHREAD].position = 0;
   sharedScaffolds[MYTHREAD].size = outputSize;
   if (sharedScaffolds[MYTHREAD].ptr == NULL) { DIE("Could not upc_alloc %ld bytes for sharedScaffolds\n", outputSize); }
   //sharedScaffoldsPointers[MYTHREAD] = (shared[] char*) upc_alloc(outputSize * sizeof(char));
   //sharedScaffoldsPositions[MYTHREAD] = 0;
   
   upc_barrier;
   
#ifdef DEBUG
   printf("Thread %d: Done with sharedScaffolds.  outputSize %ld\n", MYTHREAD, outputSize);
#endif

   /* Cache locally the directory of shared locations for scaffold report files */
    
   int q;
   ScaffoldReportFiles *cachedSharedScaffolds = NULL;
   if (MYTHREAD == 0) {
     cachedSharedScaffolds = (ScaffoldReportFiles*) malloc_chk(THREADS * sizeof(ScaffoldReportFiles));
     for (q = 0; q < THREADS; q++) {
        cachedSharedScaffolds[q] = sharedScaffolds[q];
     }
   }
   
   upc_barrier;
   
   serS = UPC_TICKS_NOW();

#ifdef DEBUG
   printf("Thread %d: Done with parallel\n", MYTHREAD);
#endif

   int64_t nGaps = 0;
   /* From now on the execution should be serialized and only THREAD 0 executes operations, process objects according to their sorted length */
   splintTracking_endTies_t *splintTracking_endTies;
   splintDataTie_t *splintTiesDataHeap;
   int64_t statsBestTie[END_TYPES];

   if (MYTHREAD == 0) {
      
      serS = UPC_TICKS_NOW();
      
      for (i=0; i < END_TYPES; i++) {
         statsBestTie[i] = 0;
      }
      
      localTieHeap = (tie_t*) malloc_chk(totalLinks * sizeof(tie_t));
      int64_t toGet, posNow = 0;
      for (i=0; i<THREADS; i++) {
         toGet = tieSharedHeapPtr[i];
         if (toGet != 0) {
            upc_memget( &localTieHeap[posNow], tieSharedHeap[i], toGet * sizeof(tie_t));
            posNow += toGet;
         }
      }
      
      totalTiesData = posNow;
      
      /* On position 2*i goes tie for Object_i.3 and on position 2*i+1 goes tie for Object_i.5 */
      endTies = (endTies_t*) malloc_chk(2 * nTotalObjects * sizeof(endTies_t));
      splintTracking_endTies = (splintTracking_endTies_t*) malloc_chk(2 * nTotalObjects * sizeof(splintTracking_endTies_t));

      tiesDataHeap = (dataTie_t*) malloc_chk(2 * totalTiesData * sizeof(dataTie_t));
      splintTiesDataHeap = (splintDataTie_t*) malloc_chk(2 * totalTiesData * sizeof(splintDataTie_t));

      for (i = 0; i < 2 * nTotalObjects; i++) {
         endTies[i].nTies = 0;
         splintTracking_endTies[i].nTies = 0;
      }
      

#ifdef DEBUG
      printf("Thread %d: initialized endTies\n", MYTHREAD);
#endif
      /* First count the entries */
      //localTieHeap = (tie_t*) remoteTieHeap;
      
      
      
      for (i=0; i<totalTiesData; i++) {
         //posInEndTies = (localTieHeap[i].ending1 == 3) ? 2*localTieHeap[i].end1_id : 2*localTieHeap[i].end1_id+1 ;
          if (i >= totalLinks)
              DIE("Index exceeds total links: %ld >= %ld\n"
                  "Check parameters: is -L set correctly?\n", i, totalLinks);
         if ( localTieHeap[i].ending1 == 3 ) {
            posInEndTies = 2*localTieHeap[i].end1_id;
         } else {
            posInEndTies = 2*localTieHeap[i].end1_id + 1;
         }
         endTies[posInEndTies].nTies++;
         splintTracking_endTies[posInEndTies].nTies++;

         //posInEndTies = (localTieHeap[i].ending2 == 3) ? 2*localTieHeap[i].end2_id : 2*localTieHeap[i].end2_id+1 ;
         if ( localTieHeap[i].ending2 == 3 ) {
            posInEndTies = 2*localTieHeap[i].end2_id;
         } else {
            posInEndTies = 2*localTieHeap[i].end2_id + 1;
         }
         endTies[posInEndTies].nTies++;
         splintTracking_endTies[posInEndTies].nTies++;

      }
      
#ifdef DEBUG
      printf("Thread %d: Done with counting endTies\n", MYTHREAD);
#endif

      /* Now do the allocation so that the data lies serially in each entry of endTies */
      runningPtr = 0;
      endTies[0].ties = &tiesDataHeap[runningPtr];
      endTies[0].cur_pos = 0;
      
      splintTracking_endTies[0].ties = &splintTiesDataHeap[runningPtr];
      splintTracking_endTies[0].cur_pos = 0;
      
      runningPtr += endTies[0].nTies;
      
      for (i=1; i<2*nTotalObjects; i++) {
         endTies[i].ties = &tiesDataHeap[runningPtr];
         
         splintTracking_endTies[i].ties = &splintTiesDataHeap[runningPtr];

         endTies[i].cur_pos = 0;
         
         splintTracking_endTies[i].cur_pos = 0;

         runningPtr += endTies[i].nTies;
      }
   
      /* Privatization of remoteLengthArray */
      localONOObjectLengths = (oNoObject*) malloc_chk(nTotalObjects * sizeof(oNoObject));
      sortedByLen = (oNoObject*) malloc_chk(nTotalObjects * sizeof(oNoObject));
      
      //memcpy(localONOObjectLengths, (oNoObject*) remoteLengthArray, nTotalObjects * sizeof(oNoObject));
      
      int64_t heapNumber, cur2, curEnt2;
      curEnt2 = (nTotalObjects+N_ACTIVE_THREADS-1)/N_ACTIVE_THREADS;
      oNoObject *tmpBuf2 = (oNoObject*) malloc_chk(curEnt2 * sizeof(oNoObject));
      
      for (heapNumber = 0; heapNumber < N_ACTIVE_THREADS; heapNumber++) {
         upc_memget(tmpBuf2, sharingScratchpad[heapNumber], curEnt2 * sizeof(oNoObject));
         cur2 = 0;
         for (i=heapNumber; i < nTotalObjects; i+= N_ACTIVE_THREADS) {
            localONOObjectLengths[i] = tmpBuf2[cur2];
            sortedByLen[i] = localONOObjectLengths[i];
            cur2++;
         }
      }
      
      free(tmpBuf2);
      
      //for (i = 0; i < THREADS; i++) {
      //   cachedSharingScratchpad[i] = sharingScratchpad[i];
      //}
      
      
      
      //for (i = 0 ; i < nTotalObjects; i++) {
      //   heapNumber = i % N_ACTIVE_THREADS;
      //   localONOObjectLengths[i] = *(cachedSharingScratchpad[heapNumber] + i/N_ACTIVE_THREADS);
      //   sortedByLen[i] = localONOObjectLengths[i];
      //}

#ifdef DEBUG
      printf("Thread %d: Ready to iterate\n", MYTHREAD);
#endif
      
      /* Now iterate over the data and store ties appropriately */
      for (i=0; i<totalTiesData; i++) {
         nGapL = localTieHeap[i].nGapLinks;
         gapE = localTieHeap[i].gapEstimate;
         gapU = localTieHeap[i].gapUncertainty;
         
         if ( localTieHeap[i].ending1 == 3 ) {
            posInEndTies = 2*localTieHeap[i].end1_id;
         } else {
            posInEndTies = 2*localTieHeap[i].end1_id+1;
         }

         cPos = endTies[posInEndTies].cur_pos;
         
         dataArray = endTies[posInEndTies].ties;
         
         splintDataArray = splintTracking_endTies[posInEndTies].ties;
         
         dataArray[cPos].nGapLinks = nGapL;
         dataArray[cPos].gapEstimate = gapE;
         dataArray[cPos].gapUncertainty = gapU;
         dataArray[cPos].end2_id = localTieHeap[i].end2_id;
         dataArray[cPos].ending2 = localTieHeap[i].ending2;
         
         splintDataArray[cPos].end2_id = localTieHeap[i].end2_id;
         splintDataArray[cPos].ending2 = localTieHeap[i].ending2;
         splintDataArray[cPos].isSplinted = localTieHeap[i].splintFlag;

         
         dataArray[cPos].endLength = localONOObjectLengths[localTieHeap[i].end2_id].length;
         endTies[posInEndTies].cur_pos++;
         splintTracking_endTies[posInEndTies].cur_pos++;
         
         
         if ( localTieHeap[i].ending2 == 3 ) {
            posInEndTies = 2*localTieHeap[i].end2_id;
         } else {
            posInEndTies = 2*localTieHeap[i].end2_id+1;
         }
         
         cPos = endTies[posInEndTies].cur_pos;
         dataArray = endTies[posInEndTies].ties;
         
         splintDataArray = splintTracking_endTies[posInEndTies].ties;

         dataArray[cPos].nGapLinks = nGapL;
         dataArray[cPos].gapEstimate = gapE;
         dataArray[cPos].gapUncertainty = gapU;
         dataArray[cPos].end2_id = localTieHeap[i].end1_id;
         dataArray[cPos].ending2 = localTieHeap[i].ending1;
         
         splintDataArray[cPos].end2_id = localTieHeap[i].end1_id;
         splintDataArray[cPos].ending2 = localTieHeap[i].ending1;
         splintDataArray[cPos].isSplinted = localTieHeap[i].splintFlag;
         
         dataArray[cPos].endLength = localONOObjectLengths[localTieHeap[i].end1_id].length;
         endTies[posInEndTies].cur_pos++;
         splintTracking_endTies[posInEndTies].cur_pos++;

      }
      
      serE = UPC_TICKS_NOW();
      storTiesTime = UPC_TICKS_TO_SECS(serE - serS);
      
#ifdef DEBUG
      printf("Thread %d: Done with iterations.. %ld\n", MYTHREAD, nTotalObjects);
#endif

      /* Now endTies array is ready to be used */
      //sortedByLen = (oNoObject*) remoteLengthArray;
      
      start_sort = UPC_TICKS_NOW();
      
      qsort(sortedByLen, nTotalObjects, sizeof(oNoObject), oNocmpfunc);
      
      end_sort = UPC_TICKS_NOW();
      sorting_time = UPC_TICKS_TO_SECS(end_sort - start_sort);

#ifdef DEBUG
      printf("Thread %d: Done with sort %ld\n", MYTHREAD, nTotalObjects);
#endif
      
      serS = UPC_TICKS_NOW();
      
      /* Thread 0 serially marks ends */
      /* On position 2*i goes tie for Object_i.3 and on position 2*i+1 goes tie for Object_i.5 */
      endMarks = (char*) malloc_chk(2 * nTotalObjects * sizeof(char));
      
      /* Marking ends */
      for (i=0; i<nTotalObjects; i++) {
         
#ifdef DEBUG
         printf("Thread %d: Marking end %ld\n", MYTHREAD, i);
#endif
      
         cur_end = sortedByLen[i].my_id;
         cur_ending = 5;
         endMark = markEnd(cur_end, cur_ending, endTies, splintTracking_endTies, localONOObjectLengths, localContigDepth, peakDepth, depthInfoAvailable, merSize, srfFlag, scaffReportLocal);
         endMarks[2*cur_end+1] = endMark;
         
         cur_ending = 3;
         endMark = markEnd(cur_end, cur_ending, endTies, splintTracking_endTies, localONOObjectLengths, localContigDepth, peakDepth, depthInfoAvailable, merSize, srfFlag, scaffReportLocal);
         endMarks[2*cur_end] = endMark;
      }

#ifdef DEBUG
      printf("Thread %d: Marked ends\n", MYTHREAD);
#endif
      
      serE = UPC_TICKS_NOW();
      markendTime = UPC_TICKS_TO_SECS(serE - serS);
      
      serS = UPC_TICKS_NOW();

      
      /* Finding best ties */
      char *suspended;
      int piece;
      char bestTie;
      int closestRes;
      char closestResEnding;
      int bestTiePos, bestTiedArrayPos;
      end_t *bestTieArray;
      endList_t *bestTiedByArray;
      endList_t *newEndObject;
      
      suspended = (char*) malloc_chk(nTotalObjects * sizeof(char));
      for (i = 0; i < nTotalObjects; i++) {
         suspended[i] = 0;
      }
      
      bestTieArray = (end_t*) malloc_chk(2 * nTotalObjects * sizeof(end_t));
      bestTiedByArray = (endList_t*) malloc_chk(2 * nTotalObjects * sizeof(endList_t));
      for (i = 0; i < 2*nTotalObjects; i++) {
         bestTieArray[i].end_id = UNDEFINED;
         bestTiedByArray[i].end_id = UNDEFINED;
         bestTiedByArray[i].next = NULL;
      }
      
#ifdef DEBUG
      printf("Thread %d: Starting totalObjects: %ld\n", MYTHREAD, nTotalObjects);
#endif

      for (i=0; i<nTotalObjects; i++) {
         piece = sortedByLen[i].my_id;
         if (piece >= nTotalObjects) 
             DIE("Too many suspended pieces. piece=%d (expecting %ld)\n", piece, nTotalObjects);
         if (suspended[piece] == 0) {
            cur_ending = 5;
            endMark = endMarks[2*piece+1];
            bestTie = endMark;
            
            statsBestTie[endMark]++;
            
            if ( endMark == UNMARKED ) {
               bestTie = bestTieFunction(piece, cur_ending, endTies, localONOObjectLengths, suspendable, endMarks, suspended, &closestRes, &closestResEnding);
               bestTiePos = 2*piece+1;
               bestTieArray[bestTiePos].end_id = closestRes;
               bestTieArray[bestTiePos].ending = closestResEnding;
               
               if (bestTie != NO_GOOD_TIES) {

                  if ( closestResEnding == 3 ) {
                     bestTiedArrayPos = 2*closestRes ;
                  } else {
                     bestTiedArrayPos = 2*closestRes + 1 ;
                  }
                  
                  if (bestTiedByArray[bestTiedArrayPos].end_id == UNDEFINED) {
                     bestTiedByArray[bestTiedArrayPos].end_id = piece;
                     bestTiedByArray[bestTiedArrayPos].ending = cur_ending;
                  } else {
                     newEndObject = (endList_t*) malloc_chk(sizeof(endList_t));
                     newEndObject->next = bestTiedByArray[bestTiedArrayPos].next;
                     bestTiedByArray[bestTiedArrayPos].next = newEndObject;
                     newEndObject->end_id = piece;
                     newEndObject->ending = cur_ending;
                  }
               }
            }
            
            cur_ending = 3;
            endMark = endMarks[2*piece];
            bestTie = endMark;
            
            statsBestTie[endMark]++;
            
            if ( endMark == UNMARKED ) {
               bestTie = bestTieFunction(piece, cur_ending, endTies, localONOObjectLengths, suspendable, endMarks, suspended, &closestRes, &closestResEnding);
               bestTiePos = 2*piece;
               bestTieArray[bestTiePos].end_id = closestRes;
               bestTieArray[bestTiePos].ending = closestResEnding;
               
               
               if (bestTie != NO_GOOD_TIES) {
                  
                  if ( closestResEnding == 3 ) {
                     bestTiedArrayPos = 2*closestRes;
                  } else {
                     bestTiedArrayPos = 2*closestRes+1;
                  }
                  
                  if (bestTiedByArray[bestTiedArrayPos].end_id == UNDEFINED) {
                     bestTiedByArray[bestTiedArrayPos].end_id = piece;
                     bestTiedByArray[bestTiedArrayPos].ending = cur_ending;
                  } else {
                     newEndObject = (endList_t*) malloc_chk(sizeof(endList_t));
                     newEndObject->next = bestTiedByArray[bestTiedArrayPos].next;
                     bestTiedByArray[bestTiedArrayPos].next = newEndObject;
                     newEndObject->end_id = piece;
                     newEndObject->ending = cur_ending;
                  }
               }
            }
         } else {
#ifdef DEBUG
            //fprintf(tiesFD, "BESTTIE: (suspended)\tContig%d\t%d\n", piece, (int) suspended[piece]);
#endif
         }
      }

#ifdef DEBUG
      printf("Thread %d: Placing suspended ties\n", MYTHREAD);
#endif
      
      serE = UPC_TICKS_NOW();
      bestTiesTime = UPC_TICKS_TO_SECS(serE - serS);
      
      serS = UPC_TICKS_NOW();
      
      /* Place suspended where possible */
      dataTie_t *ties5, *ties3;
      end_t check_fast[TIE_FAST_SIZE];
      end_t *check_slow;
      end_t *check;
      int checkFlag;
      checkFlag = 1;
      char endMark3, endMark5;
      int nTies3, nTies5;
      int64_t j,k;
      int search_end5, search_pos5;
      char search_ending5;
      end_t bt5, t5, t3, tie5, tie3;
      int exists, nValidSuspensions, aux_pos;
      
      for ( i = 0; i < nTotalObjects; i++ ) {
         piece = sortedByLen[i].my_id;
         
         if (suspended[piece] == 0) {
            continue;
         }
         
         endMark5 = endMarks[2*piece+1];
         endMark3 = endMarks[2*piece];
         
         if ( (endMark5 != UNMARKED) || (endMark3 != UNMARKED) ) {
            continue;
         }
         
         nTies5 = endTies[2*piece+1].nTies;
         if ( nTies5 > 0 ) {
            ties5 = endTies[2*piece+1].ties;
         } else {
            continue;
         }
         
         nTies3 = endTies[2*piece].nTies;
         if ( nTies3 > 0 ) {
            ties3 = endTies[2*piece].ties;
         } else {
            continue;
         }
         
         /* Decide to use preallocated array or allocate a new one */
         if (nTies5 <= TIE_FAST_SIZE) {
            checkFlag = 1;
            check = &(check_fast[0]);
         } else {
            checkFlag = 0;
            check_slow = (end_t*) malloc_chk(nTies5 * sizeof(end_t));
            check = check_slow;
         }
         
         for (j=0 ; j < nTies5; j++) {
            search_end5 = ties5[j].end2_id;
            search_ending5 = ties5[j].ending2;
            if ( search_ending5 == 3 ) {
               search_pos5 = 2*search_end5 ;
            } else {
               search_pos5 = 2*search_end5 + 1 ;
            }
            
            bt5 = bestTieArray[search_pos5];
            if (bt5.end_id == UNDEFINED) {
               check[j].end_id = UNDEFINED;
               continue;
            }
            exists = 0;
            k = 0;
            while ( (k < nTies3) && (exists == 0) ) {
               if ( (ties3[k].end2_id == bt5.end_id) && (ties3[k].ending2 == bt5.ending) ) {
                  exists = 1;
               }
               k++;
            }
            
            if (exists == 1) {
               check[j] = bt5;
            } else {
               check[j].end_id = UNDEFINED;
            }
         }
         
         nValidSuspensions = 0;
         for (j=0 ; j < nTies5; j++) {
            if (check[j].end_id == UNDEFINED) {
               continue;
            }
            t5.end_id = ties5[j].end2_id;
            t5.ending = ties5[j].ending2;
            t3 = check[j];
            
            if ( mutualUniqueBest(t5,t3, bestTieArray, bestTiedByArray, suspended) == 1 ) {
               if (nValidSuspensions == 0) {
                  tie5 = t5;
                  tie3 = t3;
               }
               nValidSuspensions++;
            }
         }
         
         if (nValidSuspensions == 1) {
            suspended[piece] = 0;
            
            if ( tie5.ending == 3) {
               aux_pos = 2*tie5.end_id;
            } else {
               aux_pos = 2*tie5.end_id + 1;
            }
            bestTieArray[aux_pos].end_id = piece;
            bestTieArray[aux_pos].ending = 5;
            bestTiedByArray[aux_pos].end_id = piece;
            bestTiedByArray[aux_pos].ending = 5;
            bestTiedByArray[aux_pos].next = NULL;
            aux_pos = 2 * piece + 1;
            bestTieArray[aux_pos].end_id = tie5.end_id;
            bestTieArray[aux_pos].ending = tie5.ending;
            bestTiedByArray[aux_pos].end_id = tie5.end_id;
            bestTiedByArray[aux_pos].ending = tie5.ending;
            bestTiedByArray[aux_pos].next = NULL;
            
            if ( tie3.ending == 3 ) {
               aux_pos = 2*tie3.end_id;
            } else {
               aux_pos = 2*tie3.end_id + 1;
            }
            
            bestTieArray[aux_pos].end_id = piece;
            bestTieArray[aux_pos].ending = 3;
            bestTiedByArray[aux_pos].end_id = piece;
            bestTiedByArray[aux_pos].ending = 3;
            bestTiedByArray[aux_pos].next = NULL;
            aux_pos = 2 * piece;
            bestTieArray[aux_pos].end_id = tie3.end_id;
            bestTieArray[aux_pos].ending = tie3.ending;
            bestTiedByArray[aux_pos].end_id = tie3.end_id;
            bestTiedByArray[aux_pos].ending = tie3.ending;
            bestTiedByArray[aux_pos].next = NULL;
            
            nInsertedSuspensions++;
         }
         
         if (checkFlag == 0) {
            free(check);
         }
         
      }

#ifdef DEBUG
      printf("Thread %d: Locking ends together\n", MYTHREAD);
#endif
      
      serE = UPC_TICKS_NOW();
      suspensionsTime = UPC_TICKS_TO_SECS(serE - serS);
      
      serS = UPC_TICKS_NOW();
      
      /*  Lock ends together with no competing ties */
      endLock_t *endLocks = (endLock_t*) malloc_chk(2 * nTotalObjects * sizeof(endLock_t));
      char piece_ending;
      end_t packedEnd1, BestTie;
      endList_t btb;
      dataTie_t tieInfo;
      int getResult, posBestTie, nTies;
      
      for (i = 0; i < 2 * nTotalObjects; i++) {
         endLocks[i].end_id = UNDEFINED;
      }
      
      for (i = 0; i < 2 * nTotalObjects; i++) {
         btb = bestTiedByArray[i];
         if ( btb.end_id == UNDEFINED ) {
            continue;
         }
         
         if ( endLocks[i].end_id != UNDEFINED ) {
            continue;
         }
         
         piece = i / 2;
         if ( (i % 2) == 0 ) {
            piece_ending = 3;
         } else {
            piece_ending = 5;
         }
         
         if ( suspended[piece] == 1 ) {
            continue;
         }
         
         BestTie = bestTieArray[i];
         if (BestTie.end_id == UNDEFINED) {
            continue;
         }
         
         nTies = 0;
         if (btb.next == NULL) {
            nTies = 1;
         }
         
         if ( (nTies == 1) && (btb.end_id == BestTie.end_id) && (btb.ending == BestTie.ending)  ) {
            packedEnd1.end_id = piece;
            packedEnd1.ending = piece_ending;
            getResult = getTieInfo(packedEnd1, BestTie, endTies, &tieInfo);
            if (getResult == 1) {
               endLocks[i].end_id = BestTie.end_id;
               endLocks[i].ending = BestTie.ending;
               endLocks[i].gap = tieInfo.gapEstimate;
               endLocks[i].gapUncertainty = tieInfo.gapUncertainty;
               if ( BestTie.ending == 3 ) {
                  posBestTie = 2*BestTie.end_id;
               } else {
                  posBestTie = 2*BestTie.end_id + 1;
               }
               endLocks[posBestTie].end_id = piece;
               endLocks[posBestTie].ending = piece_ending;
               endLocks[posBestTie].gap = tieInfo.gapEstimate;
               endLocks[posBestTie].gapUncertainty = tieInfo.gapUncertainty;
            }
         } else {
            endLocks[i].end_id = DIVERGENCE;
         }
      }
      
#ifdef DEBUG
      printf("Thread %d: Traverse ends locks\n", MYTHREAD);
#endif

      serE = UPC_TICKS_NOW();
      lockEndsTime = UPC_TICKS_TO_SECS(serE - serS);
      
      serS = UPC_TICKS_NOW();
      
      /* Traverse end locks to build scaffolds */
      /*FILE **scaffFD;
      scaffFD = (FILE**) malloc_chk(THREADS * sizeof(FILE*));
      char scaffName[MAX_FILE_PATH];
      int fileN;
      for (fileN = 0; fileN < THREADS; fileN++) {
         sprintf(scaffName, "%s_%d", outputPrefix ,fileN);
         get_rank_path(scaffName, fileN);
         scaffFD[fileN] = fopen_chk(scaffName, "w");
      }*/
      
      serE = UPC_TICKS_NOW();
      fileOpeningTime = UPC_TICKS_TO_SECS(serE - serS);
      
      //FILE **scaffFD = (FILE**) thePointer;
      
      serS = UPC_TICKS_NOW();
      
      int scaffCoord, scaffContigIndex, scaffReportPos, nextContigLen, nReportLines, visitReverse, r, actR, originalContigName, p0, p1, originalContigLength, contigStartCoord, contigEndCoord, originalGapSize, originalGapUncertainty;
      char *visitedContigs, contigOri;
      oNoObject curObject;
      scaffoldID = 0;
      char startEnd, preState, nextEnd, prevEnd;
      visitedContigs = (char*) malloc_chk(nTotalObjects * sizeof(char));
      
      char *loopCheck = (char*) malloc_chk(nTotalObjects * sizeof(char));
      memset(loopCheck, 0, nTotalObjects * sizeof(char));
      
      int *cleanupIndices = (int*) malloc_chk(nTotalObjects * sizeof(int));
      int toCleanUp = 0;
      int z;
      
      end_t end;
      int contigLen, startContig, posAux, inScaffold, nextContig, nextGap, nextGapUncertainty, prevContig;
      endLock_t next, preStateCompact;
      char scaffReportLine[MAX_LINE_SIZE], originalContigOri;
      scaf_entry_t *reportLines;
      scaf_report_t oldReport;
      shared[] scaf_entry_t *AuxreportLines;

      int *scaff_lens = malloc_chk(nTotalObjects * sizeof(int));

      for (i=0; i<nTotalObjects; i++) {
         visitedContigs[i] = 0;
      }
      
      
      cachedEntries = (scaf_entry_t**) malloc_chk(nTotalObjects * sizeof(scaf_entry_t*));
      for (i=0; i<nTotalObjects; i++ ) {
         cachedEntries[i] = NULL;
      }
      
      for (i=0; i<nTotalObjects; i++) {
         

         curObject = localONOObjectLengths[i];
         contigLen = curObject.length;
         contig = curObject.my_id;
         if (contigLen == 0) {
            continue;
         }
         startContig = contig;
         startEnd = 5;
         
         if (visitedContigs[contig] == 1) {
            continue;
         }
         
         preState = TERMINATION;
         //memset(loopCheck, 0, nTotalObjects * sizeof(char));
         
         while1S = UPC_TICKS_NOW();
         toCleanUp = 0;

         while(1) {
            if (loopCheck[startContig] == 1) {
               break;
            }
            
            loopCheck[startContig] = 1;
            
            cleanupIndices[toCleanUp] = startContig;
            toCleanUp++;
            
            end.end_id = startContig;
            end.ending = startEnd;
            if ( startEnd == 3 ) {
               posAux = 2*startContig;
            } else {
               posAux = 2*startContig + 1;
            }
            next = endLocks[posAux];
            if (next.end_id != UNDEFINED) {
               if ((next.end_id == DIVERGENCE) || (next.end_id == CONVERGENCE)) {
                  preState = NEXT;
                  preStateCompact = next;
                  break;
               } else {
                  startContig = next.end_id;
                  if ( next.ending == 3 ) {
                     startEnd = 5;
                  } else {
                     startEnd = 3;
                  }
                  preState = EXTENSION;
               }
            } else {
               preState = TERMINATION;
               break;
            }
         }
         
#ifdef DEBUG
         printf("Thread %d: cleanup loop %ld: toCleanUp %d\n", MYTHREAD, i, toCleanUp);
#endif

         /* Cleanup loopcheck without memset */
         for (z=0; z < toCleanUp; z++) {
            loopCheck[cleanupIndices[z]] = 0;
         }
         toCleanUp = 0;
         
         while1E = UPC_TICKS_NOW();
         while1Time += UPC_TICKS_TO_SECS(while1E - while1S);
   
         //memset(loopCheck, 0, nTotalObjects * sizeof(char));
         
         inScaffold = 0;
         nextContig = startContig;
         if ( startEnd == 3 ) {
            nextEnd = 5;
         } else {
            nextEnd = 3;
         }
         nextGap = 0;
         nextGapUncertainty = 0;
         prevContig = nextContig;
         prevEnd = nextEnd;
         scaffCoord = 1;
         scaffContigIndex = 1;
         scaffReportPos = 0;
         
         while2S = UPC_TICKS_NOW();

         while (1) {
            if (loopCheck[nextContig] == 1) {
               break;
            }
            
            loopCheck[nextContig] = 1;
            
            cleanupIndices[toCleanUp] = nextContig;
            toCleanUp++;
            
            visitedContigs[nextContig] = 1;
            nextContigLen = localONOObjectLengths[nextContig].length;
            contigOri = '+';
            
            if (nextEnd == 5) {
               contigOri = '-';
            }
            
            if (inScaffold == 1) {
               //upc_synci();
               sprintf(scaffReportLine, "Scaffold%d\tGAP%d\t%d\t%d\n", scaffoldID, scaffContigIndex, nextGap, nextGapUncertainty);
               writeS = UPC_TICKS_NOW();

               //fprintf(scaffFD[scaffoldID%THREADS], "%s", scaffReportLine);
               append_scaffold(scaffoldID, scaffReportLine, cachedSharedScaffolds);
               
               writeE = UPC_TICKS_NOW();
               writeTime += UPC_TICKS_TO_SECS(writeE - writeS);
               scaffCoord += nextGap;
               scaffContigIndex++;
               nGaps++;
            }
            
            if (srfFlag == 1) {
               oldReport = scaffReportLocal[nextContig];
               nReportLines = oldReport.nEntries;
               visitReverse = 0;
               if (contigOri == '-') {
                  visitReverse = 1;
               }
               
               //reportLines = (scaf_entry_t*) malloc_chk(nReportLines * sizeof(scaf_entry_t));
               //upc_memget(reportLines, oldReport.entries, nReportLines * sizeof(scaf_entry_t));
               
               if ( cachedEntries[nextContig] == NULL ) {
                  reportLines = (scaf_entry_t*) malloc_chk(nReportLines * sizeof(scaf_entry_t));
                  upc_memget(reportLines, oldReport.entries, nReportLines * sizeof(scaf_entry_t));
                  cachedEntries[nextContig] = reportLines;
               } else {
                  reportLines = cachedEntries[nextContig];
               }
               
               //AuxreportLines = oldReport.entries;
               //reportLines = (scaf_entry_t*) AuxreportLines;
               
               for ( r = 0; r < nReportLines; r++ ) {
                  //actR = (visitReverse == 0) ? r : nReportLines - r - 1;
                  if ( visitReverse == 0 ) {
                     actR = r;
                  } else {
                     actR = nReportLines - r - 1;
                  }
                  
                  if (reportLines[actR].type == CONTIG) {
                     originalContigName = reportLines[actR].entry_id;
                     p0 = reportLines[actR].f1;
                     p1 = reportLines[actR].f2;
                     originalContigLength = p1-p0+1;
                     //originalContigOri = (reportLines[actR].sign == PLUS) ? '+' : '-';
                     if ( reportLines[actR].sign == PLUS ) {
                        originalContigOri = '+';
                     } else {
                        originalContigOri = '-';
                     }
                     if (contigOri == originalContigOri) {
                        originalContigOri = '+';
                     } else {
                        originalContigOri = '-';
                     }
                     contigStartCoord = scaffCoord;
                     contigEndCoord = scaffCoord + originalContigLength - 1;
                     if (depthInfoAvailable == 1) {
                        fdepth = localContigDepth[originalContigName];
                        //printf("%d: Contig%d %f\n", __LINE__, originalContigName, fdepth);
                        //upc_synci();
                        sprintf(scaffReportLine, "Scaffold%d\tCONTIG%d\t%cContig%d\t%d\t%d\t%f\n", scaffoldID, scaffContigIndex, originalContigOri, originalContigName, contigStartCoord, contigEndCoord, fdepth );
                     } else {
                        //upc_synci();
                        sprintf(scaffReportLine, "Scaffold%d\tCONTIG%d\t%cContig%d\t%d\t%d\n", scaffoldID, scaffContigIndex, originalContigOri, originalContigName, contigStartCoord, contigEndCoord );
                     }

                     scaff_lens[scaffoldID] = contigEndCoord;

                     writeS = UPC_TICKS_NOW();
                     //fprintf(scaffFD[scaffoldID%THREADS], "%s", scaffReportLine);
                     append_scaffold(scaffoldID, scaffReportLine, cachedSharedScaffolds);
                     writeE = UPC_TICKS_NOW();
                     writeTime += UPC_TICKS_TO_SECS(writeE - writeS);
                     scaffCoord += originalContigLength;
                  } else {
                     originalGapSize = reportLines[actR].f1;
                     originalGapUncertainty = reportLines[actR].f2;
                     //upc_synci();
                     sprintf(scaffReportLine, "Scaffold%d\tGAP%d\t%d\t%d\n", scaffoldID, scaffContigIndex, originalGapSize, originalGapUncertainty );
                     nGaps++;
                     writeS = UPC_TICKS_NOW();
                     //fprintf(scaffFD[scaffoldID%THREADS], "%s", scaffReportLine);
                     append_scaffold(scaffoldID, scaffReportLine, cachedSharedScaffolds);
                     
                     writeE = UPC_TICKS_NOW();
                     writeTime += UPC_TICKS_TO_SECS(writeE - writeS);

                     scaffCoord += originalGapSize;
                     scaffContigIndex++;
                  }
               }
               
               //free(reportLines);
               
            } else {
               contigStartCoord = scaffCoord;
               contigEndCoord = scaffCoord + nextContigLen -1 ;
               
               if (depthInfoAvailable == 1) {
                  fdepth = localContigDepth[nextContig];
                  //printf("%d: Contig%d %f\n", __LINE__, nextContig, fdepth);
                  //upc_synci();
                  sprintf(scaffReportLine, "Scaffold%d\tCONTIG%d\t%cContig%d\t%d\t%d\t%f\n", scaffoldID, scaffContigIndex, contigOri, nextContig, contigStartCoord, contigEndCoord, fdepth );
               } else {
                  //upc_synci();
                  sprintf(scaffReportLine, "Scaffold%d\tCONTIG%d\t%cContig%d\t%d\t%d\n", scaffoldID, scaffContigIndex, contigOri, nextContig, contigStartCoord, contigEndCoord );
               }

               scaff_lens[scaffoldID] = contigEndCoord;

               writeS = UPC_TICKS_NOW();
               //fprintf(scaffFD[scaffoldID%THREADS], "%s", scaffReportLine);
               append_scaffold(scaffoldID, scaffReportLine, cachedSharedScaffolds);
               writeE = UPC_TICKS_NOW();
               writeTime += UPC_TICKS_TO_SECS(writeE - writeS);
               scaffCoord += nextContigLen;
            }
            
            inScaffold = 1;
            
            end.end_id = nextContig;
            end.ending = nextEnd;
            //posAux = ( nextEnd == 3 ) ? 2*nextContig : 2*nextContig+1 ;
            if ( nextEnd == 3 ) {
               posAux = 2*nextContig;
            } else {
               posAux = 2*nextContig + 1;
            }
            next = endLocks[posAux];
            if (next.end_id != UNDEFINED) {
               if ((next.end_id == DIVERGENCE) || (next.end_id == CONVERGENCE)) {
                  //preState = NEXT;
                  //preStateCompact = next;
                  // 		print STDERR "$nextContig ($nextEnd) $next\n";
                  break;
               } else {
                  
                  prevContig = nextContig;
                  prevEnd = nextEnd;
                  nextContig = next.end_id;
                  //nextEnd = ( next.ending == 3 ) ? 5 : 3 ;
                  if ( next.ending == 3  ) {
                     nextEnd = 5;
                  } else {
                     nextEnd = 3;
                  }
                  nextGap = next.gap;
                  nextGapUncertainty = next.gapUncertainty;
                  // print STDERR "[$nextGap +/- $nextGapUncertainty] $nextContig ($nextEnd) ";
               }
            } else {
               // print STDERR "$nextContig ($nextEnd) TERMINATION\n";
               break;
            }
            
         }
#ifdef DEBUG
         printf("Thread %d: done with scaffold: %d\n", MYTHREAD, scaffoldID);
#endif
         scaffoldID++;
         
         /* Cleanup loopcheck without memset */
         for (z=0; z < toCleanUp; z++) {
            loopCheck[cleanupIndices[z]] = 0;
         }
         
         while2E = UPC_TICKS_NOW();
         while2Time += UPC_TICKS_TO_SECS(while2E - while2S);
      }
      
      serE = UPC_TICKS_NOW();
      traverseLocksTime = UPC_TICKS_TO_SECS(serE - serS);
      
      serS = UPC_TICKS_NOW();
      
      //for (fileN = 0; fileN < THREADS; fileN++) {
      //   fclose(scaffFD[fileN]);
      //}
      
      serE = UPC_TICKS_NOW();
      fileOpeningTime += UPC_TICKS_TO_SECS(serE - serS);
      
      n50 = calcN50(scaff_lens, scaffoldID);
      free(scaff_lens);
   }

   //upc_synci();
   upc_barrier;

#ifdef DEBUG
   printf("Thread %d: ready to write files\n", MYTHREAD);
#endif
   
   if (MYTHREAD == 0) {
      for (q = 0; q < THREADS; q++) {
         sharedScaffolds[q] = cachedSharedScaffolds[q];
      }
   }

   
   upc_barrier;
   
   // sharedScaffolds may have been re-allocated to a different thread (0)
   // so save the file from the owning thread
   for(int thread = 0; thread < THREADS; thread++) {
      if (upc_threadof(sharedScaffolds[thread].ptr) == MYTHREAD) {
         char scaffName[MAX_FILE_PATH];
         sprintf(scaffName, "%s/%s_%d", base_dir, outputPrefix, thread);
         get_rank_path(scaffName, thread);
         FILE *myOutFD = fopen_chk(scaffName, "w");
         assert( upc_threadof(sharedScaffolds[thread].ptr) == MYTHREAD);
         char *myOutput = (char*) sharedScaffolds[thread].ptr;
         int64_t mySize = sharedScaffolds[thread].position;
         if (mySize >= sharedScaffolds[thread].size) 
             DIE("output is too large %ld maximum %ld\n", mySize, sharedScaffolds[thread].size);
         myOutput[mySize] = '\0';
         fprintf(myOutFD, "%s", myOutput);
         fclose(myOutFD);
#ifdef DEBUG
         printf("Thread %d: wrote file for %d\n", MYTHREAD, thread);
#endif

      }
   }
   
   upc_barrier;
   
   end_timer = UPC_TICKS_NOW();
   double locking_ties_time = UPC_TICKS_TO_SECS(end_timer - start_timer);
   
   char filenameForNscaffolds[100];
   FILE *fdforNmetaScaffs;
   
   endt = UPC_TICKS_NOW();
   if (MYTHREAD == 0) {
      printf("\nTime for executing oNo is: %f seconds\n" ,(UPC_TICKS_TO_SECS(endt-start)));
      printf("In total we found %d scaffolds\n", scaffoldID);
      printf("Inserted %ld suspensions\n", nInsertedSuspensions);
      printf("N50 is %d\n", n50);
      
      printf("BEST TIE SUMMARY\n");
      printf("======================================\n");
      printf("Accepted :\t %ld\n", statsBestTie[UNMARKED]);
      printf("======================================\n");
      printf("Rejected SELF_TIE :\t %ld\n", statsBestTie[SELF_TIE]);
      printf("Rejected TIE_COLLISION :\t %ld\n", statsBestTie[TIE_COLLISION]);
      printf("Rejected DEPTH_DROP :\t %ld\n", statsBestTie[DEPTH_DROP]);
      printf("Rejected NO_TIES :\t %ld\n", statsBestTie[NO_TIES]);
      printf("Rejected MULTIPLE_BEST_TIE :\t %ld\n", statsBestTie[MULTIPLE_BEST_TIE]);
      printf("======================================\n");
      printf("======================================\n\n\n");

      printf("Sorting objects time is: %f seconds \n", sorting_time);
      printf("Locking ties (serial) time is: %f seconds\n", locking_ties_time);

      printf("======================================\n");
      printf("Storing ties locally time is : %f\n", storTiesTime);
      printf("Marking end time is : %f\n", markendTime);
      printf("Adding suspension time is : %f\n", suspensionsTime);
      printf("Finding best ties time is : %f\n", bestTiesTime);
      printf("Locking end ties time is : %f\n", lockEndsTime);
      printf("Traversing locks time is : %f\n", traverseLocksTime);
      printf("File handling time time is : %f\n", fileOpeningTime);
      printf("======================================\n\n");
      
      printf("While loop 1 time is %f\n", while1Time);
      printf("While loop 2 time is %f\n", while2Time);
      printf("Writing output time is %f\n", writeTime);
      
      printf("======================================\n\n");
      
      sprintf(filenameForNscaffolds, "nScaffolds_%s.txt", outputPrefix);
      get_rank_path(filenameForNscaffolds, -1);
      fdforNmetaScaffs = fopen_chk(filenameForNscaffolds, "w");
      fprintf(fdforNmetaScaffs, "%d\n", scaffoldID);
      fclose(fdforNmetaScaffs);
      
      sprintf(filenameForNscaffolds, "nGaps_%s.txt", outputPrefix);
      get_rank_path(filenameForNscaffolds, -1);
      fdforNmetaScaffs = fopen_chk(filenameForNscaffolds, "w");
      fprintf(fdforNmetaScaffs, "%ld\n", nGaps);
      fclose(fdforNmetaScaffs);

      char fnameForN50[MAX_FILE_PATH];
      sprintf(fnameForN50, "n50_%s.txt", outputPrefix);
      get_rank_path(fnameForN50, -1);
      FILE *fdForN50 = fopen_chk(fnameForN50, "w");
      fprintf(fdForN50, "%d %d\n", pairThreshold, n50);
      fclose(fdForN50);
   }
   
   upc_barrier;
   
   if (MYTHREAD == 0) {
      if (srfFlag == 1) {
         free(scaffReportLocal);
         for (i=0; i<nTotalObjects; i++) {
            if (cachedEntries[i] != NULL) {
               free(cachedEntries[i]);
            }
         }
      }
      free(localContigDepth);
      free(localONOObjectLengths);
      free(sortedByLen);
   }
   
   /* Use another thread id, eg 480 to measure more aacurately communication times */
   if (MYTHREAD == THREADS / 2) {
      printf("Reading/processing inputs time is: %f seconds\n", read_input_time);
      printf("Storing links time is: %f seconds\n", storing_links_time);
      printf("Local link hash table consolidation time is: %f seconds\n", local_linktable_time );
      printf("Build ties time is: %f seconds ( %f sec for communication )\n", build_ties_time, build_ties_communication_time );
      printf("Propagating object info time is: %f seconds\n", propagating_lengths_time);
   }
   
   upc_barrier;
   
   serial_printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
                 ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));

   return 0;
}

