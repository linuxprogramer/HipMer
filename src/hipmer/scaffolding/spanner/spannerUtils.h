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
#define UNDEFINED (-1)
#define ANCHOR 0
#define OUTGAP 1
#define INTGAP 2

#define MAX_LINE_SIZE 1000

/* Struct to store aligning info */
/* Maximum number of alignments is read_length - kmer_length + 1 */
/* For human short reads 101bp long and k=51 MAX_ALIGN = 51  */
#define MAX_ALIGN 51

double write_io_time = 0.0;

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

typedef struct align_info {
   char *query_name;
   int subjects_matched;
   int query_length;
   int aligned_query_locs[2*MAX_ALIGN];
   int aligned_subject_ids[MAX_ALIGN];
   int aligned_subject_lengths[MAX_ALIGN];
   int aligned_subject_locs[2*MAX_ALIGN];
   int aligned_strand[MAX_ALIGN];
   int startStatus[MAX_ALIGN];
   int endStatus[MAX_ALIGN];
   int location[MAX_ALIGN];
} align_info;

/* contigScaffoldMap data structure */
typedef struct contigScaffoldMap_t contigScaffoldMap_t;
struct contigScaffoldMap_t {
   int cStrand;
   int scaffID;
   int sStart;
   int sEnd;
};

int assignType(int type, char *resultStr)
{
   
   if (type == ANCHOR)
      sprintf(resultStr, "ANCHOR");
   else if (type == OUTGAP)
      sprintf(resultStr, "OUTGAP");
   else if (type == INTGAP)
      sprintf(resultStr, "INTGAP");
   
   return 0;
}

int assignStat(int stat, char *resultStr)
{
   
   if (stat == GAP)
      sprintf(resultStr, "GAP");
   else if (stat == INC)
      sprintf(resultStr, "INC");
   else if (stat == TRUNC)
      sprintf(resultStr, "TRUNC");
   else if (stat == FULL)
      sprintf(resultStr, "FULL");
   
   return 0;
}

int assignDir(int dir, char *resultStr)
{
   
   if (dir == IN)
      sprintf(resultStr, "IN");
   else if (dir == OUT)
      sprintf(resultStr, "OUT");
   else if (dir == MID)
      sprintf(resultStr, "MID");
   
   return 0;
}

/* Process a pair */
int processPair(align_info *lane1, align_info *lane2, FILE *outFD, int innieRemoval, int minEndSeparation, int srfInput, int64_t *uninf, int64_t *sing) {
   int aligns1 = lane1->subjects_matched;
   int aligns2 = lane2->subjects_matched;
   int sStart1, sStart2, sStop1, sStop2, sLength1, sLength2, s1Outie, s1Innie, s2Outie, s2Innie, innieSep, outieSep, type1, type2, contigEnd1, contigEnd2;
   int64_t i, pos1 = 0, pos2 = 0;
   int min, contig1, contig2, strand1, strand2, startStat1, startStat2, stopStat1, stopStat2, dirStat1 ,dirStat2, innieMaxSeparation, qStart1, qStop1, qStart2, qStop2, qLength1, qLength2;
   char type1str[8];
   char type2str[8];
   char typestring[17];
   char contig1str[50];
   char contig2str[50];
   char alignInfo1string[200];
   char alignInfo2string[200];
   char startStat1string[20];
   char startStat2string[20];
   char stopStat1string[20];
   char stopStat2string[20];
   char dirStat1string[20];
   char dirStat2string[20];
   char subjectType[20];
   char strandString1[7];
   char strandString2[7];
   
   UPC_TICK_T start_io, end_io;

   
   if ((aligns1 == 0) || (aligns2 == 0)) {
      /* This read is a singleton */
      (*sing)++;
      return 0;
   }
   
   /* Find the alignments with the "earlier starts" */
   min = lane1->aligned_query_locs[0];
   for (i=1; i < aligns1; i++) {
      if (lane1->aligned_query_locs[2*i] < min) {
         min = lane1->aligned_query_locs[2*i];
         pos1 = i;
      }
   }
   
   min = lane2->aligned_query_locs[0];
   for (i=1; i < aligns2; i++) {
      if (lane2->aligned_query_locs[2*i] < min) {
         min = lane2->aligned_query_locs[2*i];
         pos2 = i;
      }
   }
   
   contig1 = lane1->aligned_subject_ids[pos1];
   strand1 = lane1->aligned_strand[pos1];
   qStart1 = lane1->aligned_query_locs[2*pos1];
   qStop1 = lane1->aligned_query_locs[2*pos1+1];
   qLength1 = lane1->query_length;
   contig2 = lane2->aligned_subject_ids[pos2];
   strand2 = lane2->aligned_strand[pos2];
   qStart2 = lane2->aligned_query_locs[2*pos2];
   qStop2 = lane2->aligned_query_locs[2*pos2+1];
   qLength2 = lane2->query_length;
   sStart1 = lane1->aligned_subject_locs[2*pos1];
   sStop1 = lane1->aligned_subject_locs[2*pos1+1];
   sLength1 = lane1->aligned_subject_lengths[pos1];
   sStart2 = lane2->aligned_subject_locs[2*pos2];
   sStop2 = lane2->aligned_subject_locs[2*pos2+1];
   sLength2 = lane2->aligned_subject_lengths[pos2];
   
   if (contig1 == contig2) {
      /* Uninformative */
      (*uninf)++;
      return 0;
   }
   
   if (innieRemoval) {
      innieMaxSeparation = 800;
      s1Innie = (strand1 == PLUS) ? sStop1 : sLength1-sStart1+1;
      s2Innie = (strand2 == PLUS) ? sStop2 : sLength2-sStart2+1;
      innieSep = s1Innie + s2Innie;
      if (innieSep < innieMaxSeparation) {
         /* Potential innie */
         return 0;
      }
   }
   
   if (minEndSeparation != 0) {
      s1Outie = (strand1 == PLUS) ? sLength1-sStart1+1 : sStop1 ;
      s2Outie = (strand2 == PLUS) ? sLength2-sStart2+1 : sStop2 ;
      outieSep = s1Outie + s2Outie;
      if (outieSep < minEndSeparation) {
         /* Potential shorty */
         return 0;
      }
   }
   
   startStat1 = lane1->startStatus[pos1];
   stopStat1 = lane1->endStatus[pos1];
   dirStat1 = lane1->location[pos1];
   startStat2 = lane2->startStatus[pos2];
   stopStat2 = lane2->endStatus[pos2];
   dirStat2 = lane2->location[pos2];
   
   type1 = ANCHOR;
   if (startStat1 == GAP) {
      type1 = OUTGAP;
   } else if (stopStat1 == GAP) {
      type1 = INTGAP;
   }
   
   type2 = ANCHOR;
   if (startStat2 == GAP) {
      type2 = OUTGAP;
   } else if (stopStat2 == GAP) {
      type2 = INTGAP;
   }
   
   contigEnd1 = (strand1 == PLUS) ? 3 : 5 ;
   contigEnd2 = (strand2 == PLUS) ? 3 : 5 ;
   
   /* Print PAIR string */
   if (srfInput)
      sprintf(subjectType, "Scaffold");
   else
      sprintf(subjectType, "Contig");
   assignType(type1, type1str);
   assignType(type2, type2str);
   sprintf(typestring, "%s.%s", type1str, type2str);
   sprintf(contig1str, "%s%d.%d", subjectType, contig1, contigEnd1);
   sprintf(contig2str, "%s%d.%d", subjectType, contig2, contigEnd2);
   assignStat(startStat1, startStat1string);
   assignStat(startStat2, startStat2string);
   assignStat(stopStat1, stopStat1string);
   assignStat(stopStat2, stopStat2string);
   assignDir(dirStat1, dirStat1string);
   assignDir(dirStat2, dirStat2string);
   if (strand1 == PLUS)
      sprintf(strandString1, "Plus");
   else
      sprintf(strandString1, "Minus");
   if (strand2 == PLUS)
      sprintf(strandString2, "Plus");
   else
      sprintf(strandString2, "Minus");
   sprintf(alignInfo1string, "[%s.%s.%s %s %d %d %d %s%d %d %d %d %s]", startStat1string, stopStat1string, dirStat1string, lane1->query_name ,qStart1, qStop1, qLength1, subjectType ,contig1, sStart1, sStop1, sLength1, strandString1);
   sprintf(alignInfo2string, "[%s.%s.%s %s %d %d %d %s%d %d %d %d %s]", startStat2string, stopStat2string, dirStat2string, lane2->query_name ,qStart2, qStop2, qLength2, subjectType, contig2, sStart2, sStop2, sLength2, strandString2);

   
   start_io = UPC_TICKS_NOW();

   fprintf(outFD, "PAIR\t%s\t%s\t%s\t\%s\t%s\n", typestring, contig1str, alignInfo1string, contig2str, alignInfo2string);
   
   end_io = UPC_TICKS_NOW();
   write_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
   
   return 1;

}

int splitSrfLine(char *srfLine, int *scaffoldId, int *contigId, int *sStart, int *sEnd, int *cStrand) {
   char *token;
   char *aux;
   token = strtok_r(srfLine, "\t", &aux);
   (*scaffoldId) = atoi(token+8);
   token = strtok_r(NULL, "\t", &aux);
   if ((*token) != 'C')
      return FAIL;
   token = strtok_r(NULL, "\t", &aux);
   (*cStrand) = ( (*token) == '+' ) ? PLUS : MINUS;
   (*contigId) = atoi(token+7);
   token = strtok_r(NULL, "\t", &aux);
   (*sStart) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*sEnd) = atoi(token);

   return SUCCESS;
}

/* Parses a single alignment */
/* Returns 1 if the alignemnt is valid. Returns 0 if the alignment is a guard alignment */
int splitAlignment(char *input_map, char *read_id, int *qStart, int *qStop, int *qLength, int *subject, int *sStart ,int *sStop, int *sLength, int *strand){
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

/* Convert to scaffold coordinate system */
int convertToScaffoldCoordinates(int qStart, int qStop, int qLength, int *subject, int *strand, int *sStart, int *sStop, int *sLength, int *projectedStart, int *projectedEnd, shared[1] contigScaffoldMap_t *contigScaffoldMap, shared[1] int64_t *scaffLengths) {
   
   int scaffLen, contigScaffStart, contigScaffEnd, scaffID, contigScaffStrand, scaffStrand, scaffStart, scaffStop;
   contigScaffoldMap_t contigScaffoldMapEntry = contigScaffoldMap[(*subject)];
   int unalignedStart = qStart - 1;
   int unalignedEnd = qLength - qStop;


   if (contigScaffoldMapEntry.scaffID == UNDEFINED) {
      return FAIL;
   } else {
      contigScaffStrand = contigScaffoldMapEntry.cStrand;
      scaffID = contigScaffoldMapEntry.scaffID;
      contigScaffStart = contigScaffoldMapEntry.sStart;
      contigScaffEnd = contigScaffoldMapEntry.sEnd;
      scaffLen = scaffLengths[scaffID];
      if ( scaffLen == 0)
         return FAIL;
      scaffStrand = (contigScaffStrand == (*strand)) ? PLUS : MINUS ;
      scaffStart = (contigScaffStrand == PLUS) ? contigScaffStart + (*sStart) - 1 : contigScaffEnd - (*sStop) + 1 ;
      scaffStop = (contigScaffStrand == PLUS) ? contigScaffStart + (*sStop) - 1 : contigScaffEnd - (*sStart) + 1 ;
      
      (*subject) = scaffID;
      (*sLength) = scaffLen;
      (*sStart) = scaffStart;
      (*strand) = scaffStrand;
      (*sStop) = scaffStop;
      (*projectedStart) = ((*strand) == PLUS) ? (*sStart) - unalignedStart : (*sStop) + unalignedStart ;
      (*projectedEnd) = ((*strand) == PLUS) ? (*sStop) + unalignedEnd :  (*sStart) - unalignedEnd;
      
      return SUCCESS;
   }
}

/* Assess alignment location relative to scaffold (pointing out, pointing in, or in the middle) */
int assessLocationAlignment(int strand, int projectedStart, int projectedEnd, int sLength, int endDistance, int *locationRes, int64_t *uninf) {
   
   int location;
   if (strand  == PLUS) {
      if (projectedStart > (sLength - endDistance)) {
         location = OUT;
      } else if ( projectedEnd < endDistance ) {
         location = IN;
      } else {
         location = MID;
      }
   } else {
      if (projectedStart < endDistance) {
         location = OUT;
      } else if ( projectedEnd > (sLength - endDistance)) {
         location = IN;
      } else {
         location = MID;
      }
   }
   
   if (location != OUT) {
      (*uninf)++;
      return FAIL;
   } else {
      (*locationRes) = location;
      return SUCCESS;
   }

}

/* Assess alignment for completeness */
int assessAlignment(int strand, int qStart, int qStop, int qLength, int sStart, int sStop, int sLength, int fivePrimeWiggleRoomIn, int threePrimeWiggleRoomIn , int truncate, int *projectedStartRes, int *projectedEndRes ,int *startStatusRes, int *endStatusRes, int64_t *five, int64_t *three) {
   
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
   } else if ( (unalignedStart < fivePrimeWiggleRoom)|| (truncate == 5) ) {
      startStatus = INC;
   } else {
      (*five)++;
      return FAIL;
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
   } else if ( (unalignedEnd < threePrimeWiggleRoom) || (truncate == 3) ) {
      endStatus = INC;
   } else {
      (*three)++;
      return FAIL;
   }
   
   (*startStatusRes) = startStatus;
   (*endStatusRes) = endStatus;
   (*projectedStartRes) = projectedStart;
   (*projectedEndRes) = projectedEnd;
   
   return SUCCESS;
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
