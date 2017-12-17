#ifndef __ONO_UTILS_H
#define __ONO_UTILS_H

#include <assert.h>
#include "../../common/upc_compatibility.h"
#include "../../common/common.h"

#define MAX_LINESIZE 2000
#define MAX_HISTO_SIZE 100000
#define MAX_SCAFF_REPORT_LENGTH 20000
#define N_SCAFF_LINES 1000
#define MAX_LIBS 10
#define BS 1
#define UNDEFINED (-999999)
#define DIVERGENCE (-1)
#define CONVERGENCE (-2)
#define TERMINATION (-3)
#define NEXT (-4)
#define EXTENSION (-5)

#define CONTIG 0
#define GAP 1
#define PLUS 0
#define MINUS 1
#define MAX_TIES_LOCAL 1000

#define END_TYPES 6
#define UNMARKED 0
#define SELF_TIE 1
#define TIE_COLLISION 2
#define DEPTH_DROP 3
#define NO_TIES 4
#define MULTIPLE_BEST_TIE 5

#define SPLINT 0
#define SPAN 1

#define EXP_FACTOR 2
#define MAX_NAME_SIZE 255
#define MAX_LIB_NAME 10
#define MAX_NLIBS 50

#define NO_GOOD_TIES 0
#define CLOSEST_LARGE 1
#define CLOSEST_EXTENDABLE 2
#define CLOSEST 3

typedef shared[] double* sharedDoublePtr;
typedef shared[] char* sharedCharPtr;

int64_t __atoi64 (const char *nptr)
{
   int c;
   int64_t value;
   int sign;
   
   while (isspace((int)*nptr))
      ++nptr;
   
   c = (int)*nptr++;
   sign = c;
   if (c == '-' || c == '+')
      c = (int)*nptr++;
   
   value = 0;
   
   while (isdigit(c))
   {
      value = 10 * value + (c - '0');
      c = (int)*nptr++;
   }
   
   if (sign == '-')
      return -value;
   else
      return value;
}

typedef struct libInfoType {
   int insertSize;
   int stdDev;
   char linkFileNamePrefix[MAX_NAME_SIZE];
} libInfoType;

int splitLinkMetaLine(char *line, libInfoType *libInfo) {
   char *token;
   char *aux;
   
   token = strtok_r(line, "\t", &aux);
   libInfo->insertSize = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   libInfo->stdDev = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   int len = strlen(token);
   // make sure to drop the \n at the end
   strncpy(libInfo->linkFileNamePrefix, token, len-1);
   libInfo->linkFileNamePrefix[len-1] = '\0';
   return 1;
}

typedef struct tie_t {
   int end1_id;
   char ending1;
   int nGapLinks;
   int gapEstimate;
   double gapUncertainty;
   int end2_id;
   char ending2;
   char splintFlag;
} tie_t;

typedef struct dataTie_t {
   int nGapLinks;
   int gapEstimate;
   int gapUncertainty;
   int end2_id;
   char ending2;
   int endLength;
} dataTie_t;

typedef struct endTies_t {
   int nTies;
   int cur_pos;
   dataTie_t *ties;
} endTies_t;

typedef struct splintDataTie_t {
   int end2_id;
   char ending2;
   char isSplinted;
} splintDataTie_t;

typedef struct splintTracking_endTies_t {
   int nTies;
   int cur_pos;
   splintDataTie_t *ties;
} splintTracking_endTies_t;

typedef struct end_t {
   int end_id;
   char ending;
} end_t;

typedef struct endLock_t {
   int end_id;
   char ending;
   int gap;
   double gapUncertainty;
} endLock_t;

struct endList_t {
   int end_id;
   char ending;
   struct endList_t *next;
};

typedef struct endList_t endList_t;

typedef shared[] tie_t* sharedTiePtr;

typedef struct splints_stack_entry {
   int nUsedSplints;
   int gapEstimate;
   char linkFileID;
} splints_stack_entry;

typedef struct spans_stack_entry {
   int nUsedSpans;
   int gapEstimate;
   double gapUncertainty;
   char linkFileID;
} spans_stack_entry;

typedef struct lib_info {
   char lib_id;
   int insertSize;
   int stdDev;
   char linkFile[MAX_LINESIZE];
} lib_info;

typedef struct scaf_entry_t {
   char type;
   char sign;
   int entry_id;
   int f1;
   int f2;
   //int intDepth;
   double fDepth;
} scaf_entry_t;

typedef struct scaf_report_t {
   int realBases;
   int depth;
   double fdepth;
   int length;
   int nEntries;
   shared[] scaf_entry_t *entries;
} scaf_report_t;

typedef struct oNoObject{
   int length;
   int my_id;
} oNoObject;

typedef struct _ScaffoldReportFiles {
       sharedCharPtr ptr;
       int64_t position;
       int64_t size;
} ScaffoldReportFiles;

void append_scaffold(int64_t scaffoldID, char * scaffReportLine,  ScaffoldReportFiles *cachedSharedScaffolds) {
    assert(MYTHREAD == 0);
    int64_t curLen = strlen(scaffReportLine);
    int64_t destThread = scaffoldID % THREADS;
    int64_t appendSize = curLen * sizeof(char);
    //printf("destThread %ld: %s", destThread, scaffReportLine);
    if (cachedSharedScaffolds[destThread].position + appendSize >= cachedSharedScaffolds[destThread].size) {
        //printf("Thread %ld needs more memory: %ld >= %ld\n", 
        //       destThread, cachedSharedScaffolds[destThread].position + appendSize, cachedSharedScaffolds[destThread].size);
        // Unlucky thread needs more memory.
        int64_t newSize = (cachedSharedScaffolds[destThread].size + appendSize) * 3 / 2;
        sharedCharPtr newPtr = (sharedCharPtr) upc_alloc( newSize );
        assert(upc_threadof(newPtr) == MYTHREAD);
        if (newPtr == NULL) 
            DIE("append_scaffold could not re-allocate %ld bytes for unlucky thread %ld\n", newSize, destThread); 
        char *localPtr = (char*) newPtr;
        upc_memget(localPtr, cachedSharedScaffolds[destThread].ptr, cachedSharedScaffolds[destThread].position);
        sharedCharPtr oldPtr = cachedSharedScaffolds[destThread].ptr;
        cachedSharedScaffolds[destThread].ptr = newPtr;
        cachedSharedScaffolds[destThread].size = newSize;
        upc_free(oldPtr);
    }
    assert(cachedSharedScaffolds[destThread].position + appendSize < cachedSharedScaffolds[destThread].size);
    upc_memput( cachedSharedScaffolds[destThread].ptr + cachedSharedScaffolds[destThread].position, scaffReportLine, appendSize);
    cachedSharedScaffolds[destThread].position += curLen;
}



typedef shared[] oNoObject* sharedoNoObjectPtr;
typedef shared[] scaf_report_t* sharedScaffPtr;
typedef shared[] scaf_entry_t* sharedScaffEntriesPtr;


int oNocmpfunc (const void * a, const void * b)
{
   if ( ((oNoObject*) b)->length == ((oNoObject*) a)->length ) {
      return (((oNoObject*) b)->my_id - ((oNoObject*) a)->my_id);
   }
   
   return ( ((oNoObject*) b)->length - ((oNoObject*) a)->length );
}

/* Split a linkData entry */
int split_linkData(char *line, lib_info *libInfo, int *libs_no){
   char *token, *aux;
   int libnum = 0;
   char libs_string[MAX_LINESIZE];
   char inserts_string[MAX_LINESIZE];
   char std_string[MAX_LINESIZE];
   char linkFile[MAX_LINESIZE];

   token = strtok_r(line, "\t", &aux);
   strcpy(libs_string, token);
   token = strtok_r(NULL, "\t", &aux);
   strcpy(inserts_string, token);
   token = strtok_r(NULL, "\t", &aux);
   strcpy(std_string, token);
   token = strtok_r(NULL, "\t", &aux);
   strcpy(linkFile, token);
   //linkFile[strlen(linkFile)-3] = '\0';
   
   token = strtok_r(libs_string, ",", &aux);
   while (token != NULL) {
      libInfo[libnum].lib_id = *(token+3);
      strcpy(libInfo[libnum].linkFile, linkFile);
      libnum++;
      token = strtok_r(NULL, ",", &aux);
   }
   
   token = strtok_r(inserts_string, ",", &aux);
   libnum = 0;
   while (token != NULL) {
      libInfo[libnum].insertSize = atoi(token);
      libnum++;
      token = strtok_r(NULL, ",", &aux);
   }
   
   token = strtok_r(std_string, ",", &aux);
   libnum = 0;
   while (token != NULL) {
      libInfo[libnum].stdDev = atoi(token);
      libnum++;
      token = strtok_r(NULL, ",", &aux);
   }
   
   (*libs_no) = libnum;
   
   return 1;
}

int cmpFuncTie (const void * a, const void * b)
{
   const dataTie_t *p1 = (dataTie_t *)a;
   const dataTie_t *p2 = (dataTie_t *)b;
   
   if (p1->gapEstimate > p2->gapEstimate) {
      return 1;
   } else if (p1->gapEstimate < p2->gapEstimate) {
      return (-1);
   } else if (p1->endLength > p2->endLength) {
      return 1;
   } else if (p1->endLength < p2->endLength) {
      return (-1);
   }
   
   return 0;
}


char markEnd(int end, char ending, endTies_t *endTies,splintTracking_endTies_t *splintTracking_endTies,oNoObject *oNoObjectlengths ,double *contigDepths, double peakDepth, int depthInfoAvailable, int merSize, int srfFlag, scaf_report_t *scaffReport) {
   
   int testContig, i, nLinks_i, gapEstimate_i, gapUncertainty_i, contig_i, contig_j, contigLen_i, start_i, start_j, end_i, uncertainty_i, uncertainty_j, end_j, overlap, excessOverlap;
   //int posInEndTies = (ending == 3) ? 2*end : 2*end+1 ;
   int posInEndTies;
   if (ending == 3) {
      posInEndTies = 2*end;
   } else {
      posInEndTies = 2*end+1;
   }
   double contigDepth_i, testContigDepth, strain;
   int nTies = endTies[posInEndTies].nTies;
   
   if ( nTies == 0) {
      return NO_TIES;
   }
   
   testContig = end;
   testContigDepth = UNDEFINED;
   
   if ( depthInfoAvailable == 1 ) {
      if (srfFlag) {
         testContigDepth = scaffReport[testContig].fdepth/scaffReport[testContig].realBases;
      } else {
         testContigDepth = contigDepths[testContig];
      }
   }
   
   int nCollisions = 0;
   int maxStrain = 3;
   double maxDepthDropoff = 0.8;
   dataTie_t *ties = endTies[posInEndTies].ties;
#ifdef DEBUG
   double expr;
#endif
   
   for (i=0; i<nTies; i++) {
      nLinks_i = ties[i].nGapLinks;
      gapEstimate_i = ties[i].gapEstimate;
      gapUncertainty_i = ties[i].gapUncertainty;
      contig_i = ties[i].end2_id;
      
      if (contig_i == testContig) {
         return SELF_TIE;
      }
      
      if ( depthInfoAvailable == 1 ) {
         if (srfFlag) {
            contigDepth_i = scaffReport[contig_i].fdepth/scaffReport[contig_i].realBases;
         } else {
            contigDepth_i = contigDepths[contig_i];
         }
         
         if ( ((1.0*(testContigDepth-contigDepth_i))/(1.0*peakDepth)) > maxDepthDropoff ) {
            return DEPTH_DROP;
         }
      }
      
      //contigLen_i = oNoObjectlengths[contig_i].length;
   }
   
   qsort(ties, nTies, sizeof(dataTie_t), cmpFuncTie);
   char contigEnd_i, contigEnd_j;
   int posToCheck = 0, nTiesToCheck = 0, z, isSplinted = 0 ;
   splintDataTie_t *splintDataArray;
   
   for (i = 0; i < nTies-1; i++) {
      start_i = ties[i].gapEstimate;
      end_i = start_i + ties[i].endLength - 1;
      uncertainty_i = ties[i].gapUncertainty;
      
      start_j = ties[i+1].gapEstimate;
      end_j = start_j + ties[i+1].endLength - 1;
      uncertainty_j = ties[i+1].gapUncertainty;
      
      overlap = end_i - start_j + 1;
      if ( overlap > merSize-2 ) {
         excessOverlap = overlap - (merSize-2);
         strain = (1.0*excessOverlap) / (1.0*(uncertainty_i + uncertainty_j));
         if (strain > (1.0*maxStrain)){
            
            contig_i = ties[i].end2_id;
            contigEnd_i = ties[i].ending2;
            if (contigEnd_i == 3) {
               contigEnd_i = 5;
            } else {
               contigEnd_i = 3;
            }
            
            contig_j = ties[i+1].end2_id;
            contigEnd_j = ties[i+1].ending2;
            
            if ( contigEnd_i == 3 ) {
               posToCheck = 2*contig_i;
            } else {
               posToCheck = 2*contig_i+1;
            }
            
            nTiesToCheck = splintTracking_endTies[posToCheck].nTies;
            splintDataArray = splintTracking_endTies[posToCheck].ties;
            isSplinted = 0;
            
            for (z = 0; z < nTiesToCheck; z++) {
               if ((splintDataArray[z].end2_id == contig_j) && (splintDataArray[z].ending2 == contigEnd_j) && (splintDataArray[z].isSplinted == 1) ) {
                  isSplinted = 1;
               }
            }
            
            if (isSplinted == 0) {
               return TIE_COLLISION;
            }
         }
      }
   }
   
   return UNMARKED;
}

char bestTieFunction(int end, char ending, endTies_t *endTies, oNoObject *oNoObjectlengths, int suspendable, char *endMarks, char *suspended, int *closestLargeRes, char *closestLargeEndingRes) {
   
   //int posInEndTies = (ending == 3) ? 2*end : 2*end+1 ;
   int posInEndTies;
   if (ending == 3) {
      posInEndTies = 2*end;
   } else {
      posInEndTies = 2*end+1;
   }
   int nTies = endTies[posInEndTies].nTies;
   
   if ( nTies == 0) {
      (*closestLargeRes) = UNDEFINED;
      return NO_GOOD_TIES;
   }
   
   dataTie_t *ties = endTies[posInEndTies].ties;
   int testContig = end;
   int testContigLength = oNoObjectlengths[testContig].length;
   int largeObject = 0;
   if  ( testContigLength > suspendable ) {
      largeObject = 1;
   }
   
   int nGoodTies, nLinks_i, gapEstimate_i, gapUncertainty_i, contig_i;
   char ending_i, endMark_i;
   int posTiedEnd_i;
   int64_t i;
   nGoodTies = 0;
   dataTie_t *goodTies = (dataTie_t*) malloc_chk(nTies * sizeof(dataTie_t));
   
   for ( i = 0; i < nTies; i++ ) {
      nLinks_i = ties[i].nGapLinks;
      gapEstimate_i = ties[i].gapEstimate;
      gapUncertainty_i = ties[i].gapUncertainty;
      contig_i = ties[i].end2_id;
      ending_i = ties[i].ending2;
      
      //posTiedEnd_i = (ending_i == 3) ? 2 * contig_i  : 2 * contig_i + 1;
      if (ending_i == 3 ) {
         posTiedEnd_i = 2 * contig_i;
      } else {
         posTiedEnd_i = 2 * contig_i + 1;
      }
      endMark_i = endMarks[posTiedEnd_i];
      
      if (endMark_i == UNMARKED) {
         goodTies[nGoodTies] = ties[i];
         nGoodTies++;
      }
   }
   
   if (nGoodTies == 0) {
      free(goodTies);
      (*closestLargeRes) = UNDEFINED;
      return NO_GOOD_TIES;
   }
   
   qsort(goodTies, nGoodTies, sizeof(dataTie_t), cmpFuncTie);
   int *localSuspended = (int*) malloc_chk(nGoodTies * sizeof(int));
   int nSus, closestLarge, closestExtendable, sContig, contigLen_i, testEnd_i;
   char closestLargeEnding, closestExtendableEnding, sEnding, end_i, otherEnd_i;
   nSus = 0;
   
   /* If a large object, return the closest large object (if one exists) small objects may be suspended between large objects */
   
   if (largeObject == 1) {
      closestLarge = UNDEFINED;
      for (i = 0; i < nGoodTies; i++) {
         sContig = goodTies[i].end2_id;
         sEnding = goodTies[i].ending2;
         contigLen_i = goodTies[i].endLength;
         if (contigLen_i > suspendable) {
            closestLarge = sContig;
            closestLargeEnding = sEnding;
            break;
         } else {
            localSuspended[nSus] = sContig;
            nSus++;
         }
      }
      
      if (closestLarge != UNDEFINED) {
#ifdef DEBUG
         if (localSuspended[i] == 5306) {
            printf("Closest large is %d\n", closestLarge);
         }
#endif
         for (i = 0; i < nSus; i++) {
            suspended[localSuspended[i]] = 1;
         }
         
         (*closestLargeRes) = closestLarge;
         (*closestLargeEndingRes) = closestLargeEnding;
         free(goodTies);
         free(localSuspended);
         return CLOSEST_LARGE;
      }
   }
   
   /* Return the closest extendable object (if one exists) unextendable objects may be suspended */
   
   nSus = 0;
   closestExtendable = UNDEFINED;
   for (i = 0; i < nGoodTies; i++) {
      contig_i = goodTies[i].end2_id;
      end_i = goodTies[i].ending2;
      //otherEnd_i = (end_i == 3) ? 5 : 3;
      if (end_i == 3)  {
         otherEnd_i = 5;
      } else {
         otherEnd_i = 3;
      }
      
      //testEnd_i = (otherEnd_i == 3) ? 2*contig_i : 2*contig_i+1;
      if (otherEnd_i == 3) {
         testEnd_i = 2*contig_i;
      } else {
         testEnd_i = 2*contig_i+1;
      }
      endMark_i = endMarks[testEnd_i];
      if (endMark_i == UNMARKED) {
         closestExtendable = contig_i;
         closestExtendableEnding = end_i;
         break;
      } else {
         localSuspended[nSus] = contig_i;
         nSus++;
      }
   }
   
   if (closestExtendable != UNDEFINED) {
      for (i = 0; i < nSus; i++) {
         suspended[localSuspended[i]] = 1;
      }
      
      (*closestLargeRes) = closestExtendable;
      (*closestLargeEndingRes) = closestExtendableEnding;
      free(goodTies);
      free(localSuspended);
      return CLOSEST_EXTENDABLE;
   }
   
   /* Just return the closest object */
   
   (*closestLargeRes) = goodTies[0].end2_id;
   (*closestLargeEndingRes) = goodTies[0].ending2;
   free(goodTies);
   free(localSuspended);
   
   return CLOSEST;
}

int mutualUniqueBest(end_t end1, end_t end2, end_t *bestTieArray, endList_t *bestTiedByArray, char *suspended ) {
   
   int piece1 = end1.end_id;
   int piece2 = end2.end_id;
   
   if ( (suspended[piece1] == 1) || (suspended[piece2] == 1) ) {
      return 0;
   }
   
   //int pos1 = ( end1.ending == 3 ) ? 2*piece1 : 2*piece1 + 1;
   //int pos2 = ( end2.ending == 3 ) ? 2*piece2 : 2*piece2 + 1;
   int pos1, pos2;
   if ( end1.ending == 3 ) {
      pos1 = 2*piece1;
   } else {
      pos1 = 2*piece1 + 1;
   }
   
   if ( end2.ending == 3 ) {
      pos2 = 2*piece2;
   } else {
      pos2 = 2*piece2 + 1;
   }
   
   int exists1 = 0;
   int exists2 = 0;
   end_t bestTie1, bestTie2;
   endList_t  bestTiedBy1, bestTiedBy2;
   
   bestTie1 = bestTieArray[pos1];
   
#ifdef DEBUG
   if (piece1 == 5306 || piece2 == 5306) {
      printf("Ammeeeen %d\n", bestTieArray[2*5306+1].end_id);
   }
#endif
   
   if (bestTie1.end_id == UNDEFINED) {
      return 0;
   }
   
   bestTie2 = bestTieArray[pos2];
   if (bestTie2.end_id == UNDEFINED) {
      return 0;
   }
   
   if ( (bestTie1.end_id != end2.end_id) || (bestTie1.ending != end2.ending) ) {
      return 0;
   }
   
   if ( (bestTie2.end_id != end1.end_id) || (bestTie2.ending != end1.ending) ) {
      return 0;
   }
   
   bestTiedBy1 = bestTiedByArray[pos1];
   if (bestTiedBy1.end_id == UNDEFINED) {
      return 0;
   }
   
   bestTiedBy2 = bestTiedByArray[pos2];
   if (bestTiedBy2.end_id == UNDEFINED) {
      return 0;
   }
   
   if ( bestTiedBy1.next != NULL ) {
      return 0;
   }
   
   if ( (bestTiedBy1.end_id != end2.end_id) || (bestTiedBy1.ending != end2.ending) ) {
      return 0;
   }
   
   if ( bestTiedBy2.next != NULL ) {
      return 0;
   }
   
   if ( (bestTiedBy2.end_id != end1.end_id) || (bestTiedBy2.ending != end1.ending) ) {
      return 0;
   }
   
   return 1;

}

int getTieInfo(end_t end1, end_t end2, endTies_t *endTies, dataTie_t *result) {

   //int pos1 = (end1.ending == 3) ? 2*end1.end_id : 2*end1.end_id+1 ;
   int pos1;
   if (end1.ending == 3) {
      pos1 = 2*end1.end_id;
   } else {
      pos1 = 2*end1.end_id+1;
   }
   int nTies = endTies[pos1].nTies;
   if (nTies == 0) {
      return 0;
   }
   
   dataTie_t *ties;
   dataTie_t tie_i;
   ties = endTies[pos1].ties;
   int i;
   
   for (i=0; i < nTies; i++) {
      tie_i = ties[i];
      if ( (tie_i.end2_id == end2.end_id) && (tie_i.ending2 == end2.ending) ) {
         (*result) = tie_i;
         return 1;
      }
   }
   
   return 0;
}

#endif
