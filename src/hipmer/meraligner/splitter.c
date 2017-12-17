#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../common/common.h"
#include "fasta.c"

#define SEGMENT_LENGTH 51


int main(int argc, char* argv[]) {
   char fastaSegment[SEGMENT_LENGTH];
   fastaSegment[SEGMENT_LENGTH-1] = '\0';
   int towrite, total_written;
   char *filename = argv[1];
   int files_no = atoi(argv[2]);
   int contigs  = atoi(argv[3]);
   FILE *inFD, *curFD;
   int current_out_id = 0;
   int contigs_read = 0;
   int contigs_per_file = (int) ceil((double)contigs / (double)files_no);
   char outputfile_name[255];
   int i, open_new_file = 0, init = 0;
   
   /* Calculate filesizes in terms of number of Contigs */
   int *filesizes = (int*) malloc_chk(files_no * sizeof(int));
   for (i = 0; i< files_no; i++) {
      filesizes[i] = 0;
   }
   
   for (i = 0; i< contigs; i++) {
      filesizes[i%files_no]++;
   }
   
   char *seq;
   char *name;
   int length;
   FASTAFILE *ffp;
   ffp = OpenFASTA(filename);
   i = 0;

   while ( ReadFASTA(ffp, &seq, &name, &length) ) {
      
      contigs_read++;
      if (contigs_read > filesizes[i]) {
         i++;
         contigs_read = 1;
         open_new_file = 1;
      } else {
         open_new_file = 0;
      }
      if (open_new_file || init == 0 ) {
         if (init == 1) fclose(curFD);
         init = 1;
         sprintf(outputfile_name,"contigs_%d" ,current_out_id);
         current_out_id++;
         curFD = fopen_chk(outputfile_name, "w+");
      }
      
      fprintf(curFD, ">%s\n", name);
      /* Print contig in FASTA FORMAT */
      total_written = 0;
      while ( total_written < length ) {
         if (total_written + SEGMENT_LENGTH-1 < length) {
            towrite = SEGMENT_LENGTH-1;
         } else {
            towrite = length - total_written;
         }
         memcpy(fastaSegment, seq+total_written, towrite * sizeof(char));
         fastaSegment[towrite] = '\0';
         fprintf(curFD, "%s\n", fastaSegment);
         total_written += towrite;

      }
      
      free(name);
      free(seq);
   }
   CloseFASTA(ffp);
   
   printf("File with the last contig is %d\n", current_out_id);
   fclose(curFD);

   return 1;
}
