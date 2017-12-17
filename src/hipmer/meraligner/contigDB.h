/* Calculate the cardinality of seeds in a contig file and the number of subcontigs generated */
#define MAX_CONTIG_SIZE 900000

int64_t seedAndContigCardinality(char *FASTAfilename, int64_t *subContigCardinality, int subContigEffectiveLength)
{
   FASTAFILE *ffp;
   int64_t resultCardinality = 0;
   int64_t subcontigCardinality = 0;
   char *seq;
   char *name;
   int length;

   ffp = OpenFASTA(FASTAfilename);
   if(ffp == NULL) {
      printf("Thread %d: Could not open or empty file: %s in seedAndContigCardinality\n", MYTHREAD, FASTAfilename);
      (*subContigCardinality) = 0;
      return 0;
   }
   while ( ReadFASTA(ffp, &seq, &name, &length) ) {
      resultCardinality += length - KMER_LENGTH + 1;
      if (!subContigEffectiveLength)
          subcontigCardinality++;
      else 
          subcontigCardinality += (int64_t) ceil((double) length/(double) subContigEffectiveLength);
      free(name);
      free(seq);
   }
   CloseFASTA(ffp);
   
   (*subContigCardinality) = subcontigCardinality;
   return resultCardinality;
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

