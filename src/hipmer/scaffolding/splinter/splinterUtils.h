#include "../../common/upc_compatibility.h"
#define PLUS 0
#define MINUS 1
#define ON 1
#define OFF 0
#define SUCCESS 1
#define FAIL 0
#define FULL 0
#define GAP 1
#define INC 2
#define TRUNC 3
#define DIVERGENT 0
#define CONVERGENT 1
#define UNK 0
#define OUT 1
#define IN 2
#define MID 3

#define MAX_LINE_SIZE 1000

/* Struct to store aligning info */
/* Maximum number of alignments is read_length - kmer_length + 1 */
/* For human short reads 101bp long and k=51 MAX_ALIGN = 51  */
#define MAX_ALIGN 51

double write_io_time = 0.0;

typedef struct align_info {
   char read_name[MAX_LINE_SIZE];
   int contigs_matched;
   int read_length;
   int aligned_read_locs[2*MAX_ALIGN];
   int aligned_contig_ids[MAX_ALIGN];
   int aligned_contig_lengths[MAX_ALIGN];
   int aligned_contig_locs[2*MAX_ALIGN];
   int aligned_strand[MAX_ALIGN];
   int startStatus[MAX_ALIGN];
   int endStatus[MAX_ALIGN];
   int location[MAX_ALIGN];
} align_info;

/* Parses a single alignment */
/* Returns 1 if the alignemnt is valid. Returns 0 if the alignment is a guard alignment */

int splitAlignment(char *input_map, char *read_id, int *qStart, int *qStop, int *qLength, int *subject, int *sStart, int *sStop, int *sLength, int *strand){
   char *token;
   char *aux;
   
   token = strtok_r(input_map, "\t", &aux);
   if ( strcmp(token, "MERALIGNER-F") == 0 ) {
      token = strtok_r(NULL, "\t", &aux);
      strcpy(read_id, token);
      return 0;
   }
   token = strtok_r(NULL, "\t", &aux);
   strcpy(read_id, token);
   token = strtok_r(NULL, "\t", &aux);
   (*qStart) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*qStop) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*qLength) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*subject) = atoi(token+6);
   token = strtok_r(NULL, "\t", &aux);
   (*sStart) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*sStop) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*sLength) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*strand) = ( strcmp(token, "Plus") == 0) ? PLUS : MINUS;
   
   return 1;
}

/* Checks if the read-to-contig alignment meets some criteria and returns appropriate characterization */
int findReadStatus(int strand, int qStart, int qStop, int qLength, int sStart, int sStop, int sLength, int fivePrimeWiggleRoomIn, int threePrimeWiggleRoomIn ,int *startStatusRes, int *endStatusRes, int *locationRes, int64_t *my_uninformative, int64_t *my_truncated) {
   
   int endDistance, location;
   int fivePrimeWiggleRoom = fivePrimeWiggleRoomIn;
   int threePrimeWiggleRoom = threePrimeWiggleRoomIn;
   int missingStartBases, startStatus, unalignedEnd, projectedEnd, endStatus, missingEndBases;
   
   int unalignedStart = qStart -1;
   int projectedStart = (strand == PLUS) ? sStart - unalignedStart : sStop + unalignedStart;
   int projectedOff = 0;
   if (projectedStart < 1)
      projectedOff = 1 - projectedStart;
   else if (projectedStart > sLength)
      projectedOff = projectedStart - sLength;
   missingStartBases = unalignedStart - projectedOff;
   
   if (unalignedStart == 0) {
      startStatus = FULL;
   } else if ( (projectedOff > 0) && (missingStartBases < fivePrimeWiggleRoom) ) {
      startStatus = GAP;
   } else if ( unalignedStart < fivePrimeWiggleRoom ) {
      startStatus = INC;
   } else {
      startStatus = TRUNC;
   }
   
   unalignedEnd = qLength - qStop;
   projectedEnd = (strand == PLUS) ? sStop + unalignedEnd : sStart - unalignedEnd;
   projectedOff = 0;
   if (projectedEnd < 1) {
      projectedOff = 1 - projectedEnd;
   } else if (projectedEnd > sLength) {
      projectedOff = projectedEnd - sLength;
   }
   missingEndBases = unalignedEnd - projectedOff;
   if (unalignedEnd == 0) {
      endStatus = FULL;
   } else if ( (projectedOff > 0) && (missingEndBases < threePrimeWiggleRoom) ) {
      endStatus = GAP;
   } else if ( unalignedEnd < threePrimeWiggleRoom ) {
      endStatus = INC;
   } else {
      endStatus = TRUNC;
   }
   
   endDistance = qLength;
   
   location = UNK;
   if ( strand == PLUS ) {
      if (projectedStart > (sLength - endDistance)) {
         location = OUT;
      } else if ( projectedEnd < endDistance ) {
         location = IN;
      } else {
         location = MID;
      }
   } else {
      if ( projectedStart < endDistance ) {
         location = OUT;
      } else if ( projectedEnd > (sLength - endDistance) ) {
         location = IN;
      } else {
         location = MID;
      }
   }
   
   (*startStatusRes) = startStatus;
   (*endStatusRes) = endStatus;
   (*locationRes) = location;
   
   if ( (startStatus == TRUNC) || (endStatus == TRUNC) ) {
      (*my_truncated)++;
      return FAIL;
   }
   
   if ( !((startStatus == GAP) || (endStatus == GAP)) ) {
      (*my_uninformative)++;
      return FAIL;
   }
   
   if ( (startStatus == GAP) || (endStatus == GAP)  ) {
      return SUCCESS;
   }
   
   return FAIL;
}

int assignStrings(int startStatus, int endStatus, int location, char **startStatus_ptr, char **endStatus_ptr, char **location_ptr,  char *inc,  char *gap,  char *trunc,  char *ins,  char *outs,  char *mid,  char *full)
{
   
   if ( startStatus == FULL ) (*startStatus_ptr) = full;
   else if ( startStatus == INC ) (*startStatus_ptr) = inc;
   else if ( startStatus == GAP ) (*startStatus_ptr) = gap;
   else if ( startStatus == TRUNC ) (*startStatus_ptr) = trunc;
   
   if ( endStatus == FULL ) (*endStatus_ptr) = full;
   else if ( endStatus == INC ) (*endStatus_ptr) = inc;
   else if ( endStatus == GAP ) (*endStatus_ptr) = gap;
   else if ( endStatus == TRUNC ) (*endStatus_ptr) = trunc;
   
   if ( location == IN ) (*location_ptr) = ins;
   else if ( location == OUT ) (*location_ptr) = outs;
   else if ( location == MID ) (*location_ptr) = mid;
   
   return 0;
}

/* Processes current splint and prints appropriate info */
int processAlignment(align_info *cur_alignments, int nAligns, int64_t *nSplintsFound, FILE *outFD)
{
   int prev_pos, cur_pos, thisContig, nextContig, thisContigEnd, nextContigEnd, thisStrand, nextStrand, thisStopStat, nextStartStat;
   char *inc="INC", *gap = "GAP", *trunc="TRUNC", *ins="IN", *outs="OUT", *mid="MID", *full = "FULL";
   char *prev_start_status, *cur_start_status, *prev_loc, *cur_loc, *prev_end_status, *cur_end_status;
   char *pl = "Plus", *prev_strand, *cur_strand, *mi = "Minus";
   
   UPC_TICK_T start_io, end_io;
   
   
   /* IMPORTANT: According to merAligner, alignments should be already sorted by the qStart index */
   /*for (prev_pos = 0; prev_pos < nAligns - 2; prev_pos++) {
      if (cur_alignments->aligned_read_locs[2*prev_pos] > cur_alignments->aligned_read_locs[(2*prev_pos+1)]) {
         printf("FALSE ASSUMPTION!!!!\n");
      }
   }*/
   
   for (prev_pos = 0; prev_pos < nAligns - 1; prev_pos++) {
      cur_pos = prev_pos + 1;
      thisContig = cur_alignments->aligned_contig_ids[prev_pos];
      nextContig = cur_alignments->aligned_contig_ids[cur_pos];
   
      if ( thisContig != nextContig) {
         thisStopStat = cur_alignments->endStatus[prev_pos];
         nextStartStat = cur_alignments->startStatus[cur_pos];
   
         if ( (thisStopStat == GAP) && (nextStartStat == GAP) ) {
            thisStrand = cur_alignments->aligned_strand[prev_pos];
            nextStrand = cur_alignments->aligned_strand[cur_pos];
            thisContigEnd = (thisStrand == PLUS) ? 3 : 5 ;
            nextContigEnd = (nextStrand == PLUS) ? 5 : 3 ;
      
            assignStrings(cur_alignments->startStatus[prev_pos], cur_alignments->endStatus[prev_pos], cur_alignments->location[prev_pos], &prev_start_status, &prev_end_status, &prev_loc, inc, gap, trunc, ins, outs, mid, full);
            prev_strand = (cur_alignments->aligned_strand[prev_pos] == PLUS) ? pl : mi ;
      
            assignStrings(cur_alignments->startStatus[cur_pos], cur_alignments->endStatus[cur_pos], cur_alignments->location[cur_pos], &cur_start_status, &cur_end_status, &cur_loc, inc, gap, trunc, ins, outs, mid, full);
            cur_strand = (cur_alignments->aligned_strand[cur_pos] == PLUS) ? pl : mi ;
            
            start_io = UPC_TICKS_NOW();
      
            fprintf(outFD, "SINGLE\tSPLINT\tContig%d.%d\t[%s.%s.%s %s %d %d %d Contig%d %d %d %d %s]\tContig%d.%d\t[%s.%s.%s %s %d %d %d Contig%d %d %d %d %s]\n",thisContig, thisContigEnd, prev_start_status, prev_end_status , prev_loc , cur_alignments->read_name, cur_alignments->aligned_read_locs[2*prev_pos], cur_alignments->aligned_read_locs[2*prev_pos+1], cur_alignments->read_length, cur_alignments->aligned_contig_ids[prev_pos], cur_alignments->aligned_contig_locs[2*prev_pos], cur_alignments->aligned_contig_locs[2*prev_pos+1], cur_alignments->aligned_contig_lengths[prev_pos], prev_strand , nextContig, nextContigEnd, cur_start_status, cur_end_status, cur_loc, cur_alignments->read_name, cur_alignments->aligned_read_locs[2*cur_pos], cur_alignments->aligned_read_locs[2*cur_pos+1], cur_alignments->read_length, cur_alignments->aligned_contig_ids[cur_pos], cur_alignments->aligned_contig_locs[2*cur_pos], cur_alignments->aligned_contig_locs[2*cur_pos+1], cur_alignments->aligned_contig_lengths[cur_pos], cur_strand);
            (*nSplintsFound) += 1;
            
            end_io = UPC_TICKS_NOW();
            write_io_time += UPC_TICKS_TO_SECS(end_io - start_io);

         }
      }
   }
   
   return 0;
}

char *sgets( char * str, int num, char **input )
{
   char *next = *input;
   int numread = 0;
   int strpos = 0;
   
   while ( numread + 1 < num && *next ) {
      int isnewline = ( *next == '\n' );
      
      str[strpos] = *next++;
      strpos++;
      numread++;
      // newline terminates the line but is included
      if ( isnewline )
         break;
   }
   
   if ( numread == 0 )
      return NULL;  // "eof"
   
   // must have hit the null terminator or end of line
   str[strpos] = '\0';  // null terminate this tring
   // set up input for next call
   *input = next;
   return str;
}
