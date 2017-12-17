#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <inttypes.h>
#include <upc_tick.h>
#include <libgen.h>
#include <zlib.h>
#include "../../common/optlist.h"
#include "../../common/common.h"
#include "../../common/Buffer.h"
#include "../../common/upc_compatibility.h"
#include "../../common/kseq.h"
#include "../../common/findDMin.h"
#include "bubbleFinder.h"
#include "bubbleHash.h"

/* Define segment length for FASTA sequences */
#ifndef SEGMENT_LENGTH
#define SEGMENT_LENGTH 51
#endif

shared int64_t bubbleID;
shared int64_t diplotigID;
shared int64_t totDiplotigLength;
shared int64_t edges;
shared int64_t indTravs;
shared int64_t totConflicts;
shared int64_t bases_from_independent;
shared int64_t isotigID;
shared int64_t McontigID;
shared int64_t bubbleContigMinDepth;
shared int64_t totalFilteredIsotigs;

#define IS_DIPLOTIG 0
#define IS_ISOTIG 1

typedef struct {
	size_t contigId;
	int depth;
	int contigType;
	size_t recordLen;
} IdLengthDepth;

void recordIdLengthDepth(int64_t contigId, char * finalSequence, int depth, int contigType, Buffer isotFDBuffer) {
    size_t pos = getLengthBuffer(isotFDBuffer);
    IdLengthDepth ild, *ptr;
    ild.contigId = contigId;
    ild.depth = depth;
    ild.contigType = contigType;
    ild.recordLen = 0; 
    writeBuffer(isotFDBuffer, &ild, sizeof(IdLengthDepth));
    printfBuffer(isotFDBuffer, ">Contig_%lld\n%s\n", (long long) contigId, finalSequence);
    size_t newPos = getLengthBuffer(isotFDBuffer);
    ptr = (IdLengthDepth*) (getStartBuffer(isotFDBuffer) + pos);
    ptr->recordLen = newPos - pos;
    //fprintf(stderr, "Thread %d: wrote isotig %lld seqlen %lld recordLen %lld (%lld)\n", MYTHREAD, (long long) ptr->contigId, (long long) strlen(finalSequence), (long long) ptr->recordLen, (long long) contigId);
}

int64_t writeFilteredIsotigs(FILE *isotFD, Buffer isotFDBuffer) {
    int64_t written = 0;
    char *bufPos;
    IdLengthDepth ild;
    //fprintf(stderr, "Thread %d: reading from isotFDBuffer length: %lld\n", MYTHREAD, (long long) getLengthBuffer(isotFDBuffer));
    assert(getPosBuffer(isotFDBuffer) == 0);
    int64_t len;
    while ( (len = readBuffer(isotFDBuffer, &ild, sizeof(IdLengthDepth))) > 0) {
        assert(len == sizeof(IdLengthDepth));
        size_t outputLen = ild.recordLen - sizeof(IdLengthDepth);
        //fprintf(stderr, "Thread %d: read Contig %lld with %d depth and %lld size\n", MYTHREAD, (long long) ild.contigId, ild.depth, (long long) outputLen);
        if (ild.contigType == IS_DIPLOTIG || ild.depth >= bubbleContigMinDepth) {
            bufPos = getCurBuffer(isotFDBuffer);
            assert(bufPos[0] == '>');
            char *id = strchr(bufPos, '\n');
            if (id == NULL) { DIE("Invalid isotFDBuffer id!: %.*s", (int) ild.recordLen, bufPos); }
            assert(*id == '\n');
            int idLen = id - bufPos + 1;
            fwrite(getCurBuffer(isotFDBuffer), 1, idLen, isotFD);
            id++; // skip newline
            assert(*id == 'A' || *id == 'C' || *id == 'G' || *id == 'T');
            printFoldedSequence(isotFD, id, outputLen - idLen);
            written++;
	}
	isotFDBuffer->pos += outputLen;
    }
    return written;
}

int main(int argc, char **argv) {
    upc_tick_t start_time = upc_ticks_now();
   
   if (MYTHREAD == 0) bubbleID = 1;
   if (MYTHREAD == 0) edges = 0;
   if (MYTHREAD == 0) diplotigID = 0;
   if (MYTHREAD == 0) isotigID = 0;
   if (MYTHREAD == 0) indTravs = 0;
   if (MYTHREAD == 0) totConflicts = 0;
   if (MYTHREAD == 0) totDiplotigLength = 0;
   if (MYTHREAD == 0) bases_from_independent = 0;
   if (MYTHREAD == 0) McontigID = 0;
   if (MYTHREAD == 0) totalFilteredIsotigs = 0;

   upc_barrier;
   int64_t my_ind = 0;
   
   /* Use getopt() to read arguments */
   extern char *optarg;
   extern int optind;
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "c:d:N:M:o:B:b:H:D:");
   
   int kmerLength;
   int nContigs;
   int my_nContigs;
   int64_t curMetaConId;
   int64_t curDiplId;
   int64_t curIsoId;
   FILE *my_nC_fd, *nC_fd;

   const char *base_dir = ".";
   char *contigFileName, *newFileName, *depthFileName;

   UPC_TICK_T start, end;
   int bin_size = 2; // bin size should be hard coded according to Eugene
   int histogram_size = 1000;
   if (MYTHREAD == 0)
       bubbleContigMinDepth = 0;
   
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'c':
            contigFileName = thisOpt->argument;
            break;
         case 'o':
            newFileName = thisOpt->argument;
            break;
         case 'd':
            depthFileName = thisOpt->argument;
            break;
         case 'N':
            nContigs = atoi(thisOpt->argument);
            break;
         case 'M':
            my_nContigs = atoi(thisOpt->argument);
            break;
         case 'B':
             base_dir = thisOpt->argument;
             break;
         case 'b':
            bin_size = atoi(thisOpt->argument);
            break;
         case 'H':
            histogram_size = atoi(thisOpt->argument);
            break;
         case 'D':
            if (MYTHREAD == 0) bubbleContigMinDepth = atoi(thisOpt->argument);
            break;
         default:
            break;
      }
      
      free(thisOpt);
   }

   {
       char my_nContigsFile[MAX_FILE_PATH];
       sprintf(my_nContigsFile, "%s/my%s_%d.txt", base_dir, contigFileName, MYTHREAD);
       get_rank_path(my_nContigsFile, MYTHREAD);
       my_nC_fd = fopen_chk(my_nContigsFile, "r");
       fgets(my_nContigsFile, MAX_FILE_PATH, my_nC_fd);
       my_nContigsFile[strlen(my_nContigsFile)-1] = '\0';
       my_nContigs = atoi(my_nContigsFile);
       fclose(my_nC_fd);
   }
   
   /* FIXME: Only for testing purposes */
   //my_nContigs = 2*nContigs/THREADS;
   printf("Thread %d: NContigs %d, myNContigs %d\n", MYTHREAD, nContigs, my_nContigs);
   
   upc_barrier;
   
   char rstrand;
   
   int64_t *local_histo = (int64_t*) malloc_chk(histogram_size * sizeof(int64_t));
   int64_t *local_histoIsotigs = (int64_t*) malloc_chk(histogram_size * sizeof(int64_t));
   int64_t *local_histoDiplotigs = (int64_t*) malloc_chk(histogram_size * sizeof(int64_t));
   
   int int_part;
   
   gzFile contigFD;
   {
       char myContigFileName[MAX_FILE_PATH];
       sprintf(myContigFileName,"%s/%s_%d.fasta", base_dir, contigFileName, MYTHREAD);
       get_rank_path(myContigFileName, MYTHREAD);
       contigFD = gzopen(myContigFileName, "r");
       if (!contigFD) { DIE("Could not open %s\n", myContigFileName); }
   }

   FILE *merDepthFD;
   {
       char myDepthFileName[MAX_FILE_PATH];
       sprintf(myDepthFileName,"%s/merDepth_%s_%d.txt", base_dir, depthFileName, MYTHREAD);
       get_rank_path(myDepthFileName, MYTHREAD);
       merDepthFD = fopen_chk(myDepthFileName, "r");
   }

   FILE *ceaFD;
   {
       char myCeaFileName[MAX_FILE_PATH];
       sprintf(myCeaFileName,"%s/%s_%d.cea", base_dir, contigFileName, MYTHREAD);
       get_rank_path(myCeaFileName, MYTHREAD);
       ceaFD = fopen_chk(myCeaFileName, "r");
   }
   
   FILE *depth_out_FD;
   {
       char myOutputDepthFileName[MAX_FILE_PATH];
       sprintf(myOutputDepthFileName,"%s/merDepth_%s_%d.txt", base_dir, newFileName, MYTHREAD);
       get_rank_path(myOutputDepthFileName, MYTHREAD);
       depth_out_FD = fopen_chk(myOutputDepthFileName, "w");
   }

   FILE *isotFD, *isotFD_raw, *mappingFD;
   Buffer isotFDBuffer = initBuffer(8192*1024);
   {
       char myIsotigsFileName[MAX_FILE_PATH];
       sprintf(myIsotigsFileName,"%s/%s_%d.fasta", base_dir, newFileName, MYTHREAD);
       get_rank_path(myIsotigsFileName, MYTHREAD);
       isotFD = fopen_chk(myIsotigsFileName, "w");

       sprintf(myIsotigsFileName,"%s/unfiltered_%s_%d.fasta", base_dir, newFileName, MYTHREAD);
       get_rank_path(myIsotigsFileName, MYTHREAD);
       isotFD_raw = fopen_chk(myIsotigsFileName, "w");
   
       sprintf(myIsotigsFileName,"%s/mapping_%s_%d.txt", base_dir, newFileName, MYTHREAD);
       get_rank_path(myIsotigsFileName, MYTHREAD);
       mappingFD = fopen_chk(myIsotigsFileName, "w");
   }
   
   cea_t curCeaEntry;
   int hairLength = 2 * KMER_LENGTH - 1;
   char prevMer[KMER_LENGTH];
   char nextMer[KMER_LENGTH];
   char nullMer[KMER_LENGTH];
   char statusField;
   int contig;
   int contigLength;
   int j;
   
   for (j=0; j<KMER_LENGTH; j++) {
      nullMer[j] = '0';
   }
   
   for (j=0; j<histogram_size; j++) {
      local_histo[j] = 0;
      local_histoIsotigs[j] = 0;
      local_histoDiplotigs[j] = 0;
   }
   
   extensions_t extensions;
   
   /* First collectively allocate contigInfo, linkertigs, bubbleMap and tipInfo arrays */
   shared[1] cea_t *contigInfo;
   contigInfo = (shared[1] cea_t*) upc_all_alloc(nContigs, sizeof(cea_t));
   
   shared[1] extensions_t *linkertigs;
   linkertigs = (shared[1] extensions_t*) upc_all_alloc(nContigs, sizeof(extensions_t));
   
   /*
   for (j=MYTHREAD; j< nContigs; j++) {
      linkertigs[j].used_flag = UNUSED_EXT;
   }*/
   
   /* PRINT FASTA LOGISTICS */
   int64_t total_written = 0, cur_length;
   int64_t towrite;
   char fastaSegment[SEGMENT_LENGTH];
   fastaSegment[SEGMENT_LENGTH-1] = '\0';
   char *seqF;
   
   shared[1] int64_t *bubbleMap;
   bubbleMap = (shared[1] int64_t*) upc_all_alloc(nContigs, sizeof(int64_t));
   
   memory_heap_t memory_heap;
   hash_table_t *tipInfo = create_hash_table(nContigs, &memory_heap, my_nContigs);
   
   char tips[2*KMER_LENGTH];
   
   /* Reading contig end data from .cea files */
   upc_barrier;
   start = UPC_TICKS_NOW();
   
   Buffer ceaBuffer = initBuffer(MAX_LINE_SIZE);
   while ( fgetsBuffer(ceaBuffer, MAX_LINE_SIZE, ceaFD) != NULL  ) {
      splitCEA(getStartBuffer(ceaBuffer), &curCeaEntry);
      resetBuffer1(ceaBuffer);
      prevMer[0] = curCeaEntry.prevBase;
      memcpy( &prevMer[1], curCeaEntry.firstMer, (KMER_LENGTH-1) * sizeof(char) );
      memcpy(nextMer, &(curCeaEntry.lastMer[1]), (KMER_LENGTH-1) * sizeof(char));
      nextMer[KMER_LENGTH-1] = curCeaEntry.nextBase;
      contigLength = curCeaEntry.contigLen;
      
      if (contigLength <= hairLength) {
         if (!( (((curCeaEntry.prevCodeL == 'F') && isACGT(curCeaEntry.prevCodeR)) && ((curCeaEntry.nextCodeR == 'F') && isACGT(curCeaEntry.nextCodeL))) ||
            (   ((curCeaEntry.prevCodeR == 'F') && isACGT(curCeaEntry.prevCodeL)) &&
             ((curCeaEntry.nextCodeL == 'F') && isACGT(curCeaEntry.nextCodeR)))    )) {
            continue;
         }
      }
      
      statusField = ( contigLength > hairLength ) ? 1 : 0 ;
      contig = curCeaEntry.contigID;
      curCeaEntry.statusField = statusField;
      contigInfo[contig] = curCeaEntry;
      
      if ( ((curCeaEntry.prevCodeL == 'F') && isACGT(curCeaEntry.prevCodeR)) || ((curCeaEntry.nextCodeR == 'F') && isACGT(curCeaEntry.nextCodeL)) ) {
         
         if ( ((curCeaEntry.prevCodeL == 'F') && isACGT(curCeaEntry.prevCodeR)) ) {
            memcpy(&(extensions.data[0]), prevMer , KMER_LENGTH * sizeof(char));
         } else {
            memcpy(&(extensions.data[0]), nullMer , KMER_LENGTH * sizeof(char));
         }
         
         if ( ((curCeaEntry.nextCodeR == 'F') && isACGT(curCeaEntry.nextCodeL)) ) {
            memcpy(&(extensions.data[KMER_LENGTH]), nextMer , KMER_LENGTH * sizeof(char));
         } else {
            memcpy(&(extensions.data[KMER_LENGTH]), nullMer , KMER_LENGTH * sizeof(char));
         }
         
         extensions.used_flag = USED_EXT;
         linkertigs[contig] = extensions;
      
      } else if ( ((curCeaEntry.prevCodeR == 'F') && isACGT(curCeaEntry.prevCodeL)) &&  ((curCeaEntry.nextCodeL == 'F') && isACGT(curCeaEntry.nextCodeR)) ) {
         
         memcpy(tips, prevMer, KMER_LENGTH * sizeof(char));
         memcpy(&(tips[KMER_LENGTH]), nextMer, KMER_LENGTH * sizeof(char));
         
         add_tip(tipInfo, tips, contig, &memory_heap);
         
         bubbleMap[contig] = 0;
      }
   }
   
   fclose(ceaFD);

   upc_barrier;
   end = UPC_TICKS_NOW();
   double ceaProcTime = UPC_TICKS_TO_SECS(end - start);

   
   /* Store contig sequence for later use */
   shared[1] sequence_t *contigSequences;
   assert(ceaBuffer != NULL);

   int contigID;
   int contigLen;
   shared[] char *curSeq;
   
   upc_barrier;
   start = UPC_TICKS_NOW();

   contigSequences = (shared[1] sequence_t*) upc_all_alloc(nContigs, sizeof(sequence_t));
   
   kseq_t *ks = kseq_init(contigFD);
   while ( kseq_read(ks) >= 0) {
      
      /* Read a contig and its length */
      contigID = atoi(ks->name.s + 7);

      contigLen = ks->seq.l;
      
      /* Stored contigs in shared memory */
      curSeq = (shared[] char*) upc_alloc(contigLen * sizeof(char));
      contigSequences[contigID].data = curSeq;
      contigSequences[contigID].seqLength = contigLen;
      memcpy( (char*) curSeq, ks->seq.s, contigLen * sizeof(char));
   }
   
   kseq_destroy(ks);
   gzclose(contigFD);

   upc_barrier;
   end = UPC_TICKS_NOW();
   double contigProcTime = UPC_TICKS_TO_SECS(end - start);

   /* Store contig mean depth for later use */
   
   int curNmers, curContigId;
   double curDepth;
   
   upc_barrier;
   start = UPC_TICKS_NOW();

   Buffer depthBuffer = ceaBuffer; // re-use ceaBuffer
   resetBuffer(depthBuffer);
   
   while ( fgetsBuffer(depthBuffer, MAX_LINE_SIZE, merDepthFD) != NULL  ) {
      splitDEPTH(getStartBuffer(depthBuffer), &curContigId, &curNmers, &curDepth);
      contigInfo[curContigId].nMers = curNmers;
      contigInfo[curContigId].meanDepth = curDepth;
      resetBuffer(depthBuffer);
   }
   
   fclose(merDepthFD);
   freeBuffer(depthBuffer);

   upc_barrier;
   end = UPC_TICKS_NOW();
   double merdepthProcTime = UPC_TICKS_TO_SECS(end - start);
   
   /* Find bubbles from convergent tips */
   shared[] list_t* cur_tip_info;
   shared[] data_t* cur_data;
   data_t cur_data_struct;
   list_t local_tip_info;
   int64_t curBubbleId;
   
   shared[1] bubbletig_t *bubbletigs;
   bubbletigs = (shared[1] bubbletig_t*) upc_all_alloc((nContigs+1), sizeof(bubbletig_t));
   int64_t i;
   int cur_data_count;
   char rcFlag;
   
   upc_barrier;
   start = UPC_TICKS_NOW();
   
#ifdef DEBUG
   char First[KMER_LENGTH+1];
   char Last[KMER_LENGTH+1];
   
   First[KMER_LENGTH] = '\0';
   Last[KMER_LENGTH] = '\0';
#endif
   
   for (i = MYTHREAD; i < nContigs; i += THREADS) {
      cur_tip_info = tipInfo->table[i].head;
      
      while ( cur_tip_info != NULL ) {
         local_tip_info = *cur_tip_info;
         
         if (local_tip_info.nContigs > 1) {
            /* A bubble */
            curBubbleId = UPC_ATOMIC_FADD_I64(&bubbleID, 1);
            upc_memput( bubbletigs[curBubbleId].tipKey, local_tip_info.tipKey, 2 * KMER_LENGTH * sizeof(char) );
            
            cur_data = local_tip_info.data;
            cur_data_count = 0;
            
#ifdef DEBUG
            memcpy(First, local_tip_info.tipKey, KMER_LENGTH * sizeof(char));
            memcpy(Last, (local_tip_info.tipKey)+KMER_LENGTH, KMER_LENGTH * sizeof(char));
            printf("BUBBLE %lld %s.%s [", (long long) curBubbleId, First, Last);
#endif
            
            while ( cur_data != NULL) {
               cur_data_count++;
               cur_data_struct = *cur_data;
               contig = cur_data_struct.contig;
               rcFlag = cur_data_struct.rcFlag;
#ifdef DEBUG
               printf("%cContig%d ",rcFlag, contig);
#endif
               bubbleMap[contig] = curBubbleId;
               cur_data = cur_data_struct.next;
            }
            
#ifdef DEBUG
            printf("]\n");
#endif
            
            if (cur_data_count != local_tip_info.nContigs) printf("FATAL ERROR: Didn't find %d nContigs\n", local_tip_info.nContigs );
         }
         cur_tip_info = local_tip_info.next;
      }
   }
   
   upc_barrier;
   end = UPC_TICKS_NOW();
   double bubbleMapTime = UPC_TICKS_TO_SECS(end - start);
   
   /* Build the bubble-contig graph (pointsTo contains the edge info) */
   pmemory_heap_t memory_heap2;
   phash_table_t *pointsTo = create_phash_table(2*nContigs, &memory_heap2, 2*my_nContigs);
   char links[2*KMER_LENGTH];
   char rcLinks[2*KMER_LENGTH];
   char in[KMER_LENGTH];
   char rcin[KMER_LENGTH];
   char type;
   int64_t loopstart;
   loopstart = ( MYTHREAD == 0 ) ? THREADS : MYTHREAD;
   
   upc_barrier;
   start = UPC_TICKS_NOW();
   
   for (i=loopstart ; i < bubbleID; i+= THREADS ) {
      upc_memget( links, bubbletigs[i].tipKey, 2 * KMER_LENGTH * sizeof(char) );
      reverseComplementSeq(links, rcLinks, 2 * KMER_LENGTH);
      memcpy( in, links, KMER_LENGTH * sizeof(char));
      memcpy( rcin, rcLinks, KMER_LENGTH * sizeof(char));
      type = B_PLUS;
      add_edge(pointsTo, in, i, type, &memory_heap2);
      type = B_MINUS;
      add_edge(pointsTo, rcin, i, type, &memory_heap2);
      
#ifdef DEBUG
      UPC_ATOMIC_FADD_I64(&edges, 2);
#endif

   }
   
#ifdef DEBUG
   printf("Added edges from bubbletigs: %lld\n", (long long) edges);
#endif
   
   for (i = MYTHREAD; i < nContigs; i+= THREADS ) {
      if (linkertigs[i].used_flag != USED_EXT) continue;
      
      upc_memget( links, linkertigs[i].data, 2 * KMER_LENGTH * sizeof(char) );
      reverseComplementSeq(links, rcLinks, 2 * KMER_LENGTH);
      memcpy( in, links, KMER_LENGTH * sizeof(char));
      memcpy( rcin, rcLinks, KMER_LENGTH * sizeof(char));
      if (in[0] != '0') {
         type = C_PLUS;
         add_edge(pointsTo, in, i, type, &memory_heap2);
#ifdef DEBUG
         UPC_ATOMIC_FADD_I64(&edges, 1);
#endif
      }
      
      if (rcin[0] != '0') {
         type = C_MINUS;
         add_edge(pointsTo, rcin, i, type, &memory_heap2);
#ifdef DEBUG
         UPC_ATOMIC_FADD_I64(&edges, 1);
#endif
      }
   
   }
   
   upc_barrier;
   end = UPC_TICKS_NOW();
   double buildGraphTime = UPC_TICKS_TO_SECS(end - start);
   
#ifdef DEBUG
   /*Print the pointsTo data structure for debugging purposes */
   int ent;
   int oid;
   char ty;
   char ch;
   char sign;
   int y;
   char dmpPath[MAX_FILE_PATH];
   sprintf(dmpPath, "bubbleFinder.dmp");
   FILE *dmp = fopen_rank_path(dmpPath, "w", MYTHREAD);
   
   shared[] plist_t *result;
   shared[] pdata_t *iter;

   char keyb[KMER_LENGTH+1];
   keyb[KMER_LENGTH] = '\0';
   for (i=0; i < 2*nContigs; i++ ) {
      result = pointsTo->table[i].head;
      while (result != NULL){
         ent = result->nContigs;
         upc_memget(keyb, result->tipKey, KMER_LENGTH*sizeof(char));
         iter = result->data;
         oid = iter->object;
         ty = iter->type;
         if (ty == C_PLUS) { ch = 'C'; sign = '+'; };
         if (ty == C_MINUS) { ch = 'C'; sign = '-'; };
         if (ty == B_MINUS) { ch = 'B'; sign = '-'; };
         if (ty == B_PLUS) { ch = 'B'; sign = '+'; };
         
         fprintf(dmp,"%s\t%c%c%d", keyb,ch,sign,oid);
         iter = iter->next;
         while (iter != NULL) {
            oid = iter->object;
            ty = iter->type;
            fprintf(dmp,",%c%c%d",ch,sign,oid);
            iter = iter->next;
         }
         fprintf(dmp,"\n");
         result = result->next;
      }
   }
   fclose(dmp);
#endif
   
   /* Allocate array for TODO items when aborting contigs */
   
   start = UPC_TICKS_NOW();
   
   shared[1] int64_t *todo = (shared[1] int64_t*) upc_all_alloc(nContigs, sizeof(int64_t));
   
   list_t tiEntry;
   shared[] list_t *tiPtr;
   enhanced_data_t *btigs;
   int64_t nBtigs, nMers;
   int64_t g;
   shared[] data_t *data_t_ptr;
   char s;
   int64_t cid;
   double depth, *btigDepths, maxDepth;
   sequence_t cur_seq;
   char *seq;
   int cur_len, minLen, minClipLen, maxDepthIndex, btigIndex, btigLen, clipLen, diff, length_to_copy, cur_path_point_added, cur_entry_len;
   local_sequence_t *btigSeqs;
   
   /*char finalSequence[MAX_CONTIG_SIZE];
   char tailSeq[MAX_CONTIG_SIZE];
   char addBack[MAX_CONTIG_SIZE];
   char paths[MAX_CONTIG_SIZE];
   char clipSeq[MAX_CONTIG_SIZE];*/
   
   char *finalSequence, *tailSeq, *addBack, *paths, *clipSeq;
   
   finalSequence = (char*) malloc_chk(MAX_CONTIG_SIZE * sizeof(char));
   tailSeq = (char*) malloc_chk(MAX_CONTIG_SIZE * sizeof(char));
   addBack = (char*) malloc_chk(MAX_CONTIG_SIZE * sizeof(char));
   paths = (char*) malloc_chk(MAX_CONTIG_SIZE * sizeof(char));
   clipSeq = (char*) malloc_chk(MAX_CONTIG_SIZE * sizeof(char));
   
   /* Traverse the bubble-contig graph, starting from each contig */
   
   /* Allocated shared visited arrays */
   shared[1] int64_t *visitedC = (shared[1] int64_t*) upc_all_alloc(nContigs, sizeof(int64_t));
   
   for (i=MYTHREAD; i<nContigs; i+= THREADS) {
      visitedC[i] = UNVISITED;
   }
   
   shared[1] int64_t *visitedB = (shared[1] int64_t*) upc_all_alloc(bubbleID, sizeof(int64_t));
   
   for (i=MYTHREAD; i<bubbleID; i+= THREADS) {
      visitedB[i] = UNVISITED;
   }
   
   int64_t flag;
   int64_t *local_log_listC = (int64_t*) malloc_chk(nContigs * sizeof(int64_t));
   int64_t *local_log_listB = (int64_t*) malloc_chk(bubbleID * sizeof(int64_t));
   int64_t logC_ptr = 0;
   int64_t logB_ptr = 0;
   
   stream_t *upstream = (stream_t*) malloc_chk((nContigs+bubbleID) * sizeof(stream_t));
   stream_t *downstream = (stream_t*) malloc_chk((nContigs+bubbleID) * sizeof(stream_t));
   int64_t upstr_pos = 0;
   int64_t downstr_pos = 0;
   stream_t next;
   stream_t current;
   stream_t first, last;
   int getNextRes;
   int abort = 0;
   int64_t k, t;
   int64_t conflIndex;
   int64_t nTigs;
   char printMode;
   
   /* If found visited status in an entry, then go and unmark the current "visited" entries, will handle those cases later, serially. This parallel section finds cases that are traversals ending to both endpoints with "multiple connections" and finds only unvisited vertices */
   int64_t indTraversals = 0;
   int64_t conflicts = 0;
   
#ifdef DEBUG
   int64_t incons = 0;
#endif
   
   for (i=MYTHREAD; i < nContigs; i+= THREADS) {

      if (linkertigs[i].used_flag != USED_EXT) continue;
      
      flag = UPC_ATOMIC_CSWAP_I64(&(visitedC[i]), UNVISITED, VISITED);
      if ( flag == VISITED ) {
         continue; // In this case we dont have to flush any visited entries from log
      }
      abort = 0;
      logC_ptr = 0;
      logB_ptr = 0;
      local_log_listC[logC_ptr] = i;
      logC_ptr++;
      
      /* Initialize downstream */
      downstr_pos = 0;
      
      current.objectId = i;
      current.type = 'C';
      current.strand = '+';
      getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);
      
      while (getNextRes == 1) {
         
         if (next.type == 'C') {
            flag = UPC_ATOMIC_CSWAP_I64(&(visitedC[next.objectId]), UNVISITED, VISITED);
            if (flag == UNVISITED) {
               local_log_listC[logC_ptr] = next.objectId;
               logC_ptr++;
            }
         } else if (next.type == 'B') {
            flag = UPC_ATOMIC_CSWAP_I64(&(visitedB[next.objectId]), UNVISITED, VISITED);
            if (flag == UNVISITED) {
               local_log_listB[logB_ptr] = next.objectId;
               logB_ptr++;
            }
         }
         
         if ( flag == VISITED ) {
            /* We found a visited vertex... well... we have to cleanup what we visited so far */
            /* Cleanup the visitedC[] array based on the log */
            for (k = 0; k < logC_ptr; k++) {
               UPC_ATOMIC_CSWAP_I64(&(visitedC[local_log_listC[k]]), VISITED, UNVISITED);
            }
            logC_ptr = 0;
            
            /* Cleanup the visitedB[] array based on the log */
            for (k = 0; k < logB_ptr; k++) {
               UPC_ATOMIC_CSWAP_I64(&(visitedB[local_log_listB[k]]), VISITED, UNVISITED);
            }
            logB_ptr = 0;

            abort = 1;
            break;
         }
         
         /* Flag is UNVISITED, therefore we have claimed this vertex -- already added visited()=1 */
         /* Add to downstream */
         downstream[downstr_pos] = next;
         downstr_pos++;
         current = next;
         getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);

      }
      
      /* We should abort this traversal... We have already cleaned up */
      if (abort == 1) {
         
         /* All the traversal from this seed have been cleaned up! */
         conflIndex = UPC_ATOMIC_FADD_I64(&totConflicts, 1);
         todo[conflIndex] = i;
         conflicts++;
         continue;
      }
      
      current.objectId = i;
      current.type = 'C';
      current.strand = '-';
      getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);
      
      /* Initialize upstream */
      upstr_pos = 0;
      
      while (getNextRes == 1) {
         
         if (next.type == 'C') {
            flag = UPC_ATOMIC_CSWAP_I64(&(visitedC[next.objectId]), UNVISITED, VISITED);
            if (flag == UNVISITED) {
               local_log_listC[logC_ptr] = next.objectId;
               logC_ptr++;
            }
         } else if (next.type == 'B') {
            flag = UPC_ATOMIC_CSWAP_I64(&(visitedB[next.objectId]), UNVISITED, VISITED);
            if (flag == UNVISITED) {
               local_log_listB[logB_ptr] = next.objectId;
               logB_ptr++;
            }
         }
         
         if ( flag == VISITED ) {
            /* We found a visited vertex... well... we have to cleanup what we visited so far */
            /* Cleanup the visitedC[] array based on the log */
            for (k = 0; k < logC_ptr; k++) {
               UPC_ATOMIC_CSWAP_I64(&(visitedC[local_log_listC[k]]), VISITED, UNVISITED);
            }
            logC_ptr = 0;
            
            /* Cleanup the visitedB[] array based on the log */
            for (k = 0; k < logB_ptr; k++) {
               UPC_ATOMIC_CSWAP_I64(&(visitedB[local_log_listB[k]]), VISITED, UNVISITED);
            }
            logB_ptr = 0;
            
            abort = 1;
            break;
         }
         
         current = next;
         rstrand = '-';
         if (next.strand == '-') {
            rstrand = '+';
         }
         next.strand = rstrand;
         upstream[upstr_pos] = next;
         upstr_pos++;
         getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);
      
      }
      
      /* We should abort this traversal... We have already cleaned up */
      if (abort == 1) {
         
         /* All the traversal from this seed have been cleaned up! */
         conflIndex = UPC_ATOMIC_FADD_I64(&totConflicts, 1);
         todo[conflIndex] = i;
         
         conflicts++;
         continue;
      }
      
      /* Remove leading and trailing bubbles */
      /* downstream has push semantics so head is at position 0 */
      /* upstream has unshift semantics so head is at position (upstr_pos - 1) */
      
      if (upstr_pos > 0) {
         first = upstream[upstr_pos-1];
         if (first.type == 'B') {
            upstr_pos--;
         }
      }

      if (downstr_pos > 0) {
         last = downstream[downstr_pos-1];
         if (last.type == 'B') {
            downstr_pos--;
         }
      }
      
      nTigs = 1 + upstr_pos + downstr_pos;
      printMode = MAX_DEPTH;
      int64_t overall_ptr = 0;
      int64_t ind;
      stream_t tig;
      char strand;
      int64_t id;
      int64_t nContigMers = 0;
      double meanContigDepth = 0;
      int64_t nBubbles = 0;
      int64_t nCons = 0;
      indTraversals++;
      
      if (nTigs > 1) {
         //printf("Diplotig for contig %ld -- no aborted\n", i);
         curDiplId = UPC_ATOMIC_FADD_I64(&diplotigID, 1);
         nContigMers = 0;
         meanContigDepth = 0.0;
         nBubbles = 0;
         nCons = 0;
         nContigMers = 0;
         finalSequence[0] = '\0';
         
         /* Clever way to iterate over the array: upstream + "C+contigID" + downstream */
         for ( overall_ptr = 0; overall_ptr < nTigs; overall_ptr++ ) {
            if (overall_ptr < upstr_pos ) {
               tig = upstream[upstr_pos - overall_ptr - 1];
            }
            if (overall_ptr == upstr_pos) {
               tig.type = 'C';
               tig.strand = '+';
               tig.objectId = i;
            }
            if (overall_ptr > upstr_pos) {
               tig = downstream[overall_ptr - upstr_pos -1];
            }
            
            type = tig.type;
            strand = tig.strand;
            id = tig.objectId;
            
            if (type == 'B') {
               nBubbles++;
               upc_memget( links, bubbletigs[id].tipKey, 2*KMER_LENGTH*sizeof(char) );
               tiPtr = lookup_and_copy_tipInfo(tipInfo, links, &tiEntry);
               nBtigs = tiEntry.nContigs;
               if (nBtigs != 2) {
                  printf("Warning: Only two-path bubbles are currently supported\n");
               }
               
               btigs = (enhanced_data_t*) malloc_chk(nBtigs * sizeof(enhanced_data_t));
               //btigDepths = (double*) malloc_chk(nBtigs * sizeof(double));
               //btigSeqs = (local_sequence_t*) malloc_chk(nBtigs * sizeof(local_sequence_t));

               data_t_ptr = tiEntry.data;
               data_t tmp;
               for (g=0; g<nBtigs; g++) {
                  tmp = (*data_t_ptr);
                  btigs[g].contig = tmp.contig;
                  btigs[g].rcFlag = tmp.rcFlag;
                  data_t_ptr = data_t_ptr->next;
               }
               
               /* foreach my $btig (@btigs) */
               for (g=0; g<nBtigs; g++) {
                  s = btigs[g].rcFlag;
                  cid = btigs[g].contig;
                  contigInfo[cid].statusField = 2;
                  depth = contigInfo[cid].meanDepth;
                  //btigDepths[g] = depth;
                  btigs[g].depth = depth;
                  
                  cur_seq = contigSequences[cid];
                  cur_len = cur_seq.seqLength;
                  seq = (char*) malloc_chk((cur_len+1) * sizeof(char));
                  seq[cur_len] = '\0';
                  assert(cur_seq.data != NULL);
                  upc_memget(seq, cur_seq.data, (cur_len) * sizeof(char));
                  
                  if (s != strand) {
                     reverseComplementINPLACE(seq, cur_len);
                  }
                  
                  btigs[g].seq = seq;
                  btigs[g].length = cur_len;
                  //btigSeqs[g].data = seq;
                  //btigSeqs[g].seqLength = cur_len;
               }
              
               /* sort array using qsort functions */
               qsort(btigs, nBtigs, sizeof(enhanced_data_t), struct_cmp_by_length);
               
               minLen = btigs[0].length;
               minClipLen = minLen - (2*KMER_LENGTH - 4);
               tailSeq[0] = '\0';
               paths[0] = '\0';
               cur_path_point_added = 0;
               cur_entry_len = 0;
               
               if (minClipLen < 0) {
                  assert(strlen(finalSequence)+minClipLen > 0);
                  strcpy(tailSeq, &finalSequence[strlen(finalSequence)+minClipLen]);
                  finalSequence[strlen(finalSequence)+minClipLen] = '\0';
               }
               
               maxDepthIndex = 0;
               maxDepth = 0.0;
               btigIndex = 0;
               
               for (g=0; g < nBtigs; g++) {
   
                  seq = btigs[g].seq;
                  btigLen = btigs[g].length;
                  clipLen = btigLen - (2*KMER_LENGTH - 4);
                  if ( clipLen < 0 ) {
                     diff = clipLen - minClipLen;
                     memcpy(addBack, tailSeq, diff * sizeof(char));
                     addBack[diff] = '\0';
                     cur_path_point_added = strlen(paths);
                     cur_entry_len = strlen(addBack);
                     strcat(paths,addBack);
                  } else {
                     memcpy(clipSeq, &seq[KMER_LENGTH-2], clipLen * sizeof(char));
                     clipSeq[clipLen] = '\0';
                     cur_path_point_added = strlen(paths);
                     cur_entry_len = strlen(tailSeq) + strlen(clipSeq);
                     strcat(paths,tailSeq);
                     strcat(paths,clipSeq);
                  }
                  
                  depth = btigs[g].depth;
                  if (depth >= maxDepth) {
                     maxDepth = depth;
                     maxDepthIndex = cur_path_point_added;
                     length_to_copy = cur_entry_len;
                  }
                  
               }
               
               /* So far support only fo maxDepth */
               if (printMode == MAX_DEPTH) {
                  paths[maxDepthIndex+length_to_copy] = '\0';
                  strcat(finalSequence, &paths[maxDepthIndex]);
               }
               
               for (g=0; g < nBtigs; g++) {
                  free(btigs[g].seq);
               }

               free(btigs);
               
               
            } else if (type == 'C') {
            
               nCons++;
               contigInfo[id].statusField = 2;
               depth = contigInfo[id].meanDepth;
               nMers = contigInfo[id].nMers;
               nContigMers += nMers;
               meanContigDepth += nMers * depth;
               
               cur_seq = contigSequences[id];
               cur_len = cur_seq.seqLength;
               seq = (char*) malloc_chk((cur_len+1) * sizeof(char));
               seq[cur_len] = '\0';
               assert(cur_seq.data != NULL);
               upc_memget(seq, cur_seq.data, (cur_len) * sizeof(char));
               
               if (strand == '-') {
                  reverseComplementINPLACE(seq, cur_len);
               }
               
               strcat(finalSequence, seq);
               free(seq);
            }
         }
         
         meanContigDepth /= nContigMers;
         curMetaConId = UPC_ATOMIC_FADD_I64(&McontigID, 1);

         fprintf(mappingFD,"Diplotig%lld ==> Contig%lld\n", (long long) curDiplId, (long long) curMetaConId);
         
         // print filtered at the end
         recordIdLengthDepth(curMetaConId, finalSequence, meanContigDepth, IS_DIPLOTIG, isotFDBuffer);

         /* Print all *tigs in 1-line FASTA format */
         fprintf(isotFD_raw, ">Contig_%ld\n", curMetaConId);
         fprintf(isotFD_raw, "%s\n", finalSequence);

         cur_length = strlen(finalSequence);
         fprintf(depth_out_FD, "Contig%lld\t%lld\t%f\n", (long long) curMetaConId, (long long) cur_length-KMER_LENGTH+1, meanContigDepth);
         
         int_part = (int) meanContigDepth;
         int_part = int_part/bin_size;
         if (int_part < histogram_size) {
            local_histo[int_part]++;
            local_histoDiplotigs[int_part]++;
         } else {
            local_histo[histogram_size-1]++;
            local_histoDiplotigs[histogram_size-1]++;
         }

         
#ifdef DEBUG
         my_ind++;
         UPC_ATOMIC_FADD_I64(&totDiplotigLength, strlen(finalSequence));
         UPC_ATOMIC_FADD_I64(&bases_from_independent, strlen(finalSequence));
#endif
      }
   }
   
   upc_barrier;
   
   end = UPC_TICKS_NOW();
   double travGraphTime = UPC_TICKS_TO_SECS(end - start);
   
   start = UPC_TICKS_NOW();


#ifdef DEBUG
   int64_t skipped = 0;
   int64_t effLinkertigs = 0;
   char linkdatacur[2*KMER_LENGTH+1];
   linkdatacur[2*KMER_LENGTH] = '\0';
   FILE *linktigs = fopen_chk("linktigs", "w");
#endif

   /* Now just THREAD 0 will initiate traversals from the unvisited seeds */
   if (MYTHREAD == 0) {
      
      /* FIXME: Remove duplicates using sorting plus linear pass */
      int64_t *local_todo = (int64_t*) malloc_chk(totConflicts * sizeof(int64_t));
      for (t=0; t < totConflicts; t++) {
         local_todo[t] = todo[t];
      }
      
      int64_t aux1;
      int64_t aux2;
      int64_t aux3;
      int64_t dups = 0;
      
      for (aux1 = 0; aux1 < totConflicts; aux1++) {
         for (aux2 = aux1 + 1; aux2 < totConflicts;) {
            if ( local_todo[aux1] == local_todo[aux2] ) {
               for ( aux3 = aux2; aux3 < totConflicts - 1; aux3++ ) {
                  local_todo[aux3] = local_todo[aux3+1];
                  dups++;
               }
               totConflicts--;
            } else {
               aux2++;
            }
         }
      }
      
      //printf("Found %ld duplicates \n", dups);
      
      qsort(local_todo, totConflicts, sizeof(int64_t), int64_comp);
      
      for (t=0; t < totConflicts; t++) {
        i = local_todo[t];
        if (linkertigs[i].used_flag != USED_EXT) continue;
#ifdef DEBUG
        upc_memget(linkdatacur, linkertigs[i].data, 2* KMER_LENGTH * sizeof(char));
        fprintf(linktigs,"%lld\t%s\n",(long long) i,linkdatacur);
        effLinkertigs++;
#endif
         flag = UPC_ATOMIC_CSWAP_I64(&(visitedC[i]), UNVISITED, VISITED);
         if ( flag == VISITED ) {
#ifdef DEBUG
            printf("Skiped entry: %lld\n", (long long) i);
            skipped++;
#endif
            continue; // In this case we dont have to flush any visited entries from log
         }
        
         /* Initialize downstream */
         downstr_pos = 0;
         current.objectId = i;
         current.type = 'C';
         current.strand = '+';
         getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);
         
         while (getNextRes == 1) {
            
            if (next.type == 'C') {
               flag = UPC_ATOMIC_CSWAP_I64(&(visitedC[next.objectId]), UNVISITED, VISITED);
            } else if (next.type == 'B') {
               flag = UPC_ATOMIC_CSWAP_I64(&(visitedB[next.objectId]), UNVISITED, VISITED);
            }
            
            if ( flag == VISITED ) {
               getNextRes = 0;
            } else {
               /* Add to downstream */
               downstream[downstr_pos] = next;
               downstr_pos++;
               current = next;
               getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);
            }
         }
        
         /* Initialize upstream */
         upstr_pos = 0;
         current.objectId = i;
         current.type = 'C';
         current.strand = '-';
         getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);
         
         while (getNextRes == 1) {
            rstrand = '-';
            if (next.strand == '-') {
               rstrand = '+';
            }
            
            if (next.type == 'C') {
               flag = UPC_ATOMIC_CSWAP_I64(&(visitedC[next.objectId]), UNVISITED, VISITED);
            } else if (next.type == 'B') {
               flag = UPC_ATOMIC_CSWAP_I64(&(visitedB[next.objectId]), UNVISITED, VISITED);
            }
            
            if ( flag == VISITED ) {
               getNextRes = 0;
            } else {
               current = next;
               next.strand = rstrand;
               upstream[upstr_pos] = next;
               upstr_pos++;
               getNextRes = getNext(current, linkertigs, bubbletigs, pointsTo, &next);
            }
         }
         
         /* Remove leading and trailing bubbles */
         /* downstream has push semantics so head is at position 0 */
         /* upstream has unshift semantics so head is at position (upstr_pos - 1) */
        
        printMode = MAX_DEPTH;
        int64_t overall_ptr = 0;
        int64_t ind;
        stream_t tig;
        char strand;
        int64_t id;
        int64_t nContigMers = 0;
        double meanContigDepth = 0.0;
        int64_t nBubbles = 0;
        int64_t nCons = 0;
        
         if (upstr_pos > 0) {
            first = upstream[upstr_pos-1];
            strand = first.strand;
            id = first.objectId;
            if (first.type == 'B') {
               upstr_pos--;
            }
         }
        
         if (downstr_pos > 0) {
            last = downstream[downstr_pos-1];
            strand = last.strand;
            id = last.objectId;
            if (last.type == 'B') {
               downstr_pos--;
            }
         }
         
         nTigs = 1 + upstr_pos + downstr_pos;
        
         if (nTigs > 1) {
            curDiplId = UPC_ATOMIC_FADD_I64(&diplotigID, 1);
            nContigMers = 0;
            meanContigDepth = 0.0;
            nBubbles = 0;
            nCons = 0;
            nContigMers = 0;
            finalSequence[0] = '\0';
            
            /* Clever way to iterate over the array: upstream + "C+contigID" + downstream */
            for ( overall_ptr = 0; overall_ptr < nTigs; overall_ptr++ ) {
               if (overall_ptr < upstr_pos ) {
                  tig = upstream[upstr_pos - overall_ptr - 1];
               }
               if (overall_ptr > upstr_pos) {
                  tig = downstream[overall_ptr - upstr_pos -1];
               }
               if (overall_ptr == upstr_pos) {
                  tig.type = 'C';
                  tig.strand = '+';
                  tig.objectId = i;
               }
               
               type = tig.type;
               strand = tig.strand;
               id = tig.objectId;
               
#ifdef DEBUG
               printf("%c%c%lld->",type, strand, (long long) id);
#endif
               
               if (type == 'B') {
                  nBubbles++;
                  upc_memget( links, bubbletigs[id].tipKey, 2*KMER_LENGTH*sizeof(char) );
                  tiPtr = lookup_and_copy_tipInfo(tipInfo, links, &tiEntry);
                  nBtigs = tiEntry.nContigs;
                  if (nBtigs != 2) {
                     printf("Warning: Only two-path bubbles are currently supported\n");
                  }
                  
                  btigs = (enhanced_data_t*) malloc_chk(nBtigs * sizeof(enhanced_data_t));
                  //btigDepths = (double*) malloc_chk(nBtigs * sizeof(double));
                  //btigSeqs = (local_sequence_t*) malloc_chk(nBtigs * sizeof(local_sequence_t));
                  
                  data_t_ptr = tiEntry.data;
                  data_t tmp;
                  for (g=0; g<nBtigs; g++) {
                     tmp = (*data_t_ptr);
                     btigs[g].contig = tmp.contig;
                     btigs[g].rcFlag = tmp.rcFlag;
                     data_t_ptr = data_t_ptr->next;
                  }
                  
                  /* foreach my $btig (@btigs) */
                  for (g=0; g<nBtigs; g++) {
                     s = btigs[g].rcFlag;
                     cid = btigs[g].contig;
                     contigInfo[cid].statusField = 2;
                     depth = contigInfo[cid].meanDepth;
                     //btigDepths[g] = depth;
                     btigs[g].depth = depth;
                     
                     cur_seq = contigSequences[cid];
                     cur_len = cur_seq.seqLength;
                     seq = (char*) malloc_chk((cur_len+1) * sizeof(char));
                     seq[cur_len] = '\0';
                     assert(cur_seq.data != NULL);
                     upc_memget(seq, cur_seq.data, (cur_len) * sizeof(char));
                     
                     if (s != strand) {
                        reverseComplementINPLACE(seq, cur_len);
                     }
                     
                     btigs[g].seq = seq;
                     btigs[g].length = cur_len;
                     //btigSeqs[g].data = seq;
                     //btigSeqs[g].seqLength = cur_len;
                  }
                  
                  /* sort array using qsort functions */
                  qsort(btigs, nBtigs, sizeof(enhanced_data_t), struct_cmp_by_length);
                  
                  minLen = btigs[0].length;
                  minClipLen = minLen - (2*KMER_LENGTH - 4);
                  tailSeq[0] = '\0';
                  paths[0] = '\0';
                  cur_path_point_added = 0;
                  cur_entry_len = 0;
                  
                  if (minClipLen < 0) {
                     assert(strlen(finalSequence)+minClipLen > 0);
                     strcpy(tailSeq, &finalSequence[strlen(finalSequence)+minClipLen]);
                     finalSequence[strlen(finalSequence)+minClipLen] = '\0';
                  }
                  
                  maxDepthIndex = 0;
                  maxDepth = 0.0;
                  btigIndex = 0;
                  
                  for (g=0; g < nBtigs; g++) {
                     
                     seq = btigs[g].seq;
                     btigLen = btigs[g].length;
                     clipLen = btigLen - (2*KMER_LENGTH - 4);
                     if ( clipLen < 0 ) {
                        diff = clipLen - minClipLen;
                        memcpy(addBack, tailSeq, diff * sizeof(char));
                        addBack[diff] = '\0';
                        cur_path_point_added = strlen(paths);
                        cur_entry_len = strlen(addBack);
                        strcat(paths,addBack);
                     } else {
                        memcpy(clipSeq, &seq[KMER_LENGTH-2], clipLen * sizeof(char));
                        clipSeq[clipLen] = '\0';
                        cur_path_point_added = strlen(paths);
                        cur_entry_len = strlen(tailSeq) + strlen(clipSeq);
                        strcat(paths,tailSeq);
                        strcat(paths,clipSeq);
                     }
                     
                     depth = btigs[g].depth;
                     if (depth >= maxDepth) {
                        maxDepth = depth;
                        maxDepthIndex = cur_path_point_added;
                        length_to_copy = cur_entry_len;
                     }
                     
                  }
                  
                  /* So far support only fo maxDepth */
                  if (printMode == MAX_DEPTH) {
                     paths[maxDepthIndex+length_to_copy] = '\0';
                     strcat(finalSequence, &paths[maxDepthIndex]);
                  }
                  
                  for (g=0; g < nBtigs; g++) {
                     free(btigs[g].seq);
                  }
                  
                  free(btigs);
                  
                  
               } else if (type == 'C') {
                  
                  nCons++;
                  contigInfo[id].statusField = 2;
                  depth = contigInfo[id].meanDepth;
                  nMers = contigInfo[id].nMers;
                  nContigMers += nMers;
                  meanContigDepth += nMers * depth;
                  
                  cur_seq = contigSequences[id];
                  cur_len = cur_seq.seqLength;
                  seq = (char*) malloc_chk((cur_len+1) * sizeof(char));
                  seq[cur_len] = '\0';
                  assert(cur_seq.data != NULL);
                  upc_memget(seq, cur_seq.data, (cur_len) * sizeof(char));
                  
                  if (strand == '-') {
                     reverseComplementINPLACE(seq, cur_len);
                  }
                  
                  strcat(finalSequence, seq);
                  free(seq);
                  
               }
            }
#ifdef DEBUG
            printf("\n");
            UPC_ATOMIC_FADD_I64(&totDiplotigLength, strlen(finalSequence));
#endif
            meanContigDepth /= nContigMers;
            curMetaConId = UPC_ATOMIC_FADD_I64(&McontigID, 1);
            
            fprintf(mappingFD,"Diplotig%lld ==> Contig%lld\n", (long long) curDiplId, (long long) curMetaConId);
            
            // print filtered at the end
            recordIdLengthDepth(curMetaConId, finalSequence, meanContigDepth, IS_DIPLOTIG, isotFDBuffer);
            
            /* Print raw results rather than FASTA format */
            fprintf(isotFD_raw, ">Contig_%lld\n", (long long) curMetaConId);
            fprintf(isotFD_raw, "%s\n", finalSequence);
            
            cur_length = strlen(finalSequence);
            fprintf(depth_out_FD, "Contig%lld\t%lld\t%f\n", (long long) curMetaConId, (long long) cur_length-KMER_LENGTH+1, meanContigDepth);
            
            int_part = (int) meanContigDepth;
            int_part = int_part/bin_size;
            if (int_part < histogram_size) {
               local_histo[int_part]++;
               local_histoDiplotigs[int_part]++;
            } else {
               local_histo[histogram_size-1]++;
               local_histoDiplotigs[histogram_size-1]++;
            }
            
         }
      }
   }

   upc_barrier;
#ifdef DEBUG
   UPC_ATOMIC_FADD_I64(&indTravs, my_ind);
   upc_barrier;
#endif
   
   end = UPC_TICKS_NOW();
   double SerialtravGraphTime = UPC_TICKS_TO_SECS(end - start);
   
   
   
   start = UPC_TICKS_NOW();
   char stat;
   int curle;
   char *printSeq;
   char *bufferStatic = (char*) malloc_chk(MAX_CONTIG_SIZE * sizeof(char));
   sequence_t cs;
   
   for (i=MYTHREAD; i < nContigs; i+= THREADS) {
      
      stat = contigInfo[i].statusField;
      if (stat == 1) {
         curIsoId = UPC_ATOMIC_FADD_I64(&isotigID, 1);
         curMetaConId = UPC_ATOMIC_FADD_I64(&McontigID, 1);
         cs = contigSequences[i];
         curle = cs.seqLength;
         if (curle > MAX_CONTIG_SIZE - 1) {
            printSeq = (char*) malloc_chk((curle+1) * sizeof(char));
         } else {
            printSeq = bufferStatic;
         }
         
         upc_memget(printSeq, cs.data, curle * sizeof(char) );
         printSeq[curle] = '\0';
         
         fprintf(mappingFD,"Isotig%lld ==> Contig%lld\n", (long long) curIsoId, (long long) curMetaConId);
         
         // print filtered at the end
         recordIdLengthDepth(curMetaConId, printSeq, contigInfo[i].meanDepth, IS_ISOTIG, isotFDBuffer);

         /* Print raw results rather than FASTA format */
         fprintf(isotFD_raw, ">Contig_%lld\n", (long long) curMetaConId);
         fprintf(isotFD_raw, "%s\n", printSeq);
         
         /* Print contig in FASTA FORMAT */
         cur_length = strlen(printSeq);
         
         fprintf(depth_out_FD, "Contig%lld\t%lld\t%f\n", (long long) curMetaConId, (long long) cur_length-KMER_LENGTH+1, contigInfo[i].meanDepth);
         
         int_part = (int) contigInfo[i].meanDepth;
         int_part = int_part/bin_size;
         if (int_part < histogram_size) {
            local_histo[int_part]++;
            local_histoIsotigs[int_part]++;
         } else {
            local_histo[histogram_size-1]++;
            local_histoIsotigs[histogram_size-1]++;
         }

         if (curle > MAX_CONTIG_SIZE - 1) {
            free(printSeq);
         }
      }
   }
   
   upc_barrier;
   
   /* Now all reduce the local histograms */
   int padded_size = ((histogram_size+THREADS-1)/THREADS) * THREADS;
   shared[1] int64_t *dist_histograms = (shared[1] int64_t*) upc_all_alloc(THREADS * padded_size, sizeof(int64_t));
   
   shared[1] int64_t *shared_result = (shared[1] int64_t*) upc_all_alloc(histogram_size, sizeof(int64_t));
   
   for (j = 0; j < histogram_size; j++) {
      dist_histograms[MYTHREAD * padded_size + j] = local_histo[j];
   }
   
   upc_barrier;
   
   int64_t cur_sum = 0;
   
   for (j = MYTHREAD; j < histogram_size; j += THREADS) {
      cur_sum = 0;
      for (i=0; i<THREADS; i++) {
         cur_sum += dist_histograms[i * padded_size + j];
      }
      shared_result[j] = cur_sum;
   }
   
   upc_barrier;
   
   FILE *fd_histo, *fd_histoIsotigs, *fd_histoDiplotigs;
   
   if (MYTHREAD == 0) {
      for (j = 0; j < histogram_size; j++) {
         local_histo[j] = shared_result[j];
      }
      
      char path[MAX_FILE_PATH];
      sprintf(path, "contig_histogram.txt");
      fd_histo = fopen_rank_path(path, "w", -1);
      
      for (j = 0; j < histogram_size; j++) {
         fprintf(fd_histo, "%d\t%lld\n",j , (long long) local_histo[j]);
      }
      
      fclose(fd_histo);
      if (bubbleContigMinDepth == 0) {
          bubbleContigMinDepth = findDMin(local_histo, histogram_size, 1);
          if (bubbleContigMinDepth < 0) {
              DIE("Could not automatically determine the BUBBLE_MIN_DEPTH_CUTOFF.  Please examine the contig_histogram.txt and set the -D bubbleContigMinDepth");
              bubbleContigMinDepth = 1;
          }
          upc_fence;
      }
      printf("Using %lld as bubble_contig_min_depth to screen out isotigs with low coverage\n", (long long) bubbleContigMinDepth);
   }
   upc_barrier;
   // Now print the filtered isotigs  >= bubbleMinContigDepth + all diplotigs
   
   size_t filteredIsoTigs = writeFilteredIsotigs(isotFD, isotFDBuffer);
   UPC_ATOMIC_FADD_I64(&totalFilteredIsotigs, filteredIsoTigs);
   
   /*  Now all reduce the local histograms for diplotigs */
   
   for (j = 0; j < histogram_size; j++) {
      dist_histograms[MYTHREAD * padded_size + j] = local_histoDiplotigs[j];
   }
   
   upc_barrier;
   if (MYTHREAD == 0) printf("Wrote %lld filteredIsoTigs\n", (long long) totalFilteredIsotigs);
   
   for (j = MYTHREAD; j < histogram_size; j += THREADS) {
      cur_sum = 0;
      for (i=0; i<THREADS; i++) {
         cur_sum += dist_histograms[i * padded_size + j];
      }
      shared_result[j] = cur_sum;
   }
   
   upc_barrier;
   
   if (MYTHREAD == 0) {
      for (j = 0; j < histogram_size; j++) {
         local_histoDiplotigs[j] = shared_result[j];
         local_histoIsotigs[j] = local_histo[j] - local_histoDiplotigs[j];
      }
      
      char path[MAX_FILE_PATH];
      sprintf(path, "contig_histogram_isotigs.txt");
      fd_histoIsotigs = fopen_rank_path(path, "w", -1);
      sprintf(path, "contig_histogram_diplotigs.txt");
      fd_histoDiplotigs = fopen_rank_path(path, "w", -1);

      
      for (j = 0; j < histogram_size; j++) {
         fprintf(fd_histoIsotigs, "%d\t%lld\n",j , (long long) local_histoIsotigs[j]);
         fprintf(fd_histoDiplotigs, "%d\t%lld\n",j , (long long) local_histoDiplotigs[j]);
      }
      
      fclose(fd_histoIsotigs);
      fclose(fd_histoDiplotigs);
   }
   
   upc_barrier;
   
   end = UPC_TICKS_NOW();
   double isotigPrintingTime = UPC_TICKS_TO_SECS(end - start);


   /* Report performance breakdown */
   if (MYTHREAD == 0) {
      printf("Processing cea files %f seconds\n", ceaProcTime);
      printf("Processing contig files %f seconds\n", contigProcTime);
      printf("Processing merDepth files %f seconds\n", merdepthProcTime);
      printf("Filling bubbleMap %f seconds\n", bubbleMapTime);
      printf("Building bubble-contig graph time %f seconds\n", buildGraphTime);
      printf("Parallel traversal bubble-contig graph time %f seconds\n", travGraphTime);
      printf("Serial traversal bubble-contig graph time %f seconds\n", SerialtravGraphTime);
      printf("Isotig printing time %f seconds\n", isotigPrintingTime);
      printf("\n=========================================\n");
      printf("Total time %f seconds\n", ceaProcTime + contigProcTime + merdepthProcTime + bubbleMapTime + buildGraphTime + travGraphTime + isotigPrintingTime);
      printf("Total bubbles: %" PRId64 "\n", (int64_t)bubbleID-1);
#ifdef DEBUG
      printf("Total ind travs is %lld, total conflicts are %lld\n", (long long) indTravs, (long long) totConflicts);
      printf("Total edges: %lld\n", (long long) edges);
      printf("Total effective linkertigs: %lld\n", (long long) effLinkertigs);
      printf("Total skipped: %lld\n", (long long) skipped);
      printf("Total inconsistencies: %lld\n",(long long)  incons);
      printf("Total independent traversals contributed %lld bases\n",(long long)  bases_from_independent);
      printf("Total length of diplotigs is %lld bases\n",(long long)  totDiplotigLength);
#endif
      printf("Total diplotigs found %lld\n",(long long)  diplotigID);
      printf("Total isotigs found %lld\n", (long long) isotigID);
      printf("Total filtered isotigs + diplotigs %lld\n", (long long) totalFilteredIsotigs);
      
      char countFileName[MAX_FILE_PATH];
      sprintf(countFileName,"n%s.txt", newFileName);
      FILE *countFD = fopen_rank_path(countFileName, "w", -1);
      fprintf(countFD, "%lld\n", (long long) diplotigID+isotigID);
      fclose(countFD);

   }
   
   fclose(isotFD);
   fclose(isotFD_raw);
   fclose(mappingFD);
   fclose(depth_out_FD);

   upc_barrier;
   
	if (!MYTHREAD)
		printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
			   ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
	return 0;
}
