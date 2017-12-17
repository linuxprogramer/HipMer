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
#define UNDEFINED -1

#define HISTOGRAM_SIZE 200000
#define MAX_LINE_SIZE 1000

/* Struct to store aligning info */
/* Maximum number of alignments is read_length - kmer_length + 1 */
/* For human short reads 101bp long and k=51 MAX_ALIGN = 51  */
#define MAX_ALIGN 51

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
   char *read_name;
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

int cmpfunc (const void * a, const void * b)
{
   return ( *(int64_t*)a - *(int64_t*)b );
}

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
int testAlignStatus(int qStart, int qStop, int qLength, int sStart, int sStop, int sLength, int strand ,int fivePrimeWiggleRoom , int threePrimeWiggleRoom, int innieRemoval, int shortPair) {
   
   int missingStartBases, startStatus, unalignedEnd, projectedEnd, endStatus, missingEndBases, minBackDistance, backDistance;

   int unalignedStart = qStart -1;
   int projectedStart = (strand == PLUS) ? sStart - unalignedStart : sStop + unalignedStart;
   
   if (innieRemoval) {
      minBackDistance = shortPair;
      backDistance = projectedStart;
      if (strand == MINUS) {
         backDistance = sLength - projectedStart;
      }
      if (backDistance <= minBackDistance) {
         return FAIL;
      }
   }
   
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
   projectedEnd = (strand ==PLUS) ? sStop + unalignedEnd : sStart - unalignedEnd;
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
   
   if (((startStatus == FULL) || (startStatus == INC)) && ( (endStatus == FULL) || (endStatus == INC) ) )
      return SUCCESS;
   else
      return FAIL;
}
