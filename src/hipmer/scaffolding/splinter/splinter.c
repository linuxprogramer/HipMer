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

shared int64_t tot_splints;
shared int64_t singletons;
shared int64_t uninformatives;
shared int64_t truncated;


double read_io_time = 0.0;

#include "../../common/upc_compatibility.h"
#include "../../common/common.h"
#include "splinterUtils.h"

int main(int argc, char **argv) {
    upc_tick_t start_time = upc_ticks_now();

   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "l:m:F:T:r:B:");
   FILE *laneFD1, *outFD;
   char *libname = NULL;
   int minMatch, fivePrimeWiggleRoom = 5, threePrimeWiggleRoom = 5;
   char outputMerAligner[MAX_FILE_PATH];
   char outputSplinter[MAX_FILE_PATH];
   char *merAlignerSuffixName;
   align_info cur_alignments;
   char *lineBuffers = (char*) calloc_chk(MAX_LINE_SIZE*3, 1);
   char *cur_read = lineBuffers;
   char *prev_read = cur_read + MAX_LINE_SIZE;
   char *align = prev_read + MAX_LINE_SIZE;
   int splitRes, qStart, qStop, qLength, subject, sStart, sStop, sLength, strand, cur_pos;
   int64_t splints_found = 0;
   UPC_TICK_T start, end;
   int testResult, startStatusRes, endStatusRes, locationRes;
   int64_t my_sigleton = 0;
   int64_t my_uninformative = 0;
   int64_t my_truncated = 0;
   int64_t my_splints = 0;
   
   UPC_TICK_T start_io, end_io;

   char *base_dir = ".";
   
   cur_read[0] = '\0';
   prev_read[0] = '\0';
   
   if (MYTHREAD == 0) {
      singletons = 0;
      uninformatives = 0;
      truncated = 0;
   }
   
   char my_nContigsFile[255];
   FILE *my_nC_fd;
   int readLength = 0;
   
   /* Process the input arguments */
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'l':
            libname = thisOpt->argument;
            break;
         case 'T':
            threePrimeWiggleRoom = atoi(thisOpt->argument);
            break;
         case 'F':
            fivePrimeWiggleRoom = atoi(thisOpt->argument);
            break;
         case 'm':
            minMatch = atoi(thisOpt->argument);
            break;
         case 'r':
            readLength = atoi(thisOpt->argument);
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
       DIE("Please specify a -l libname\n");

   upc_barrier;
   start = UPC_TICKS_NOW();

   sprintf(outputSplinter, "%s/%s-splints_%d", base_dir, libname, MYTHREAD);
   get_rank_path(outputSplinter, MYTHREAD);
   outFD = fopen_chk(outputSplinter, "w+");
   for(int readNumber = 1; readNumber <= 2; readNumber++) {
       sprintf(outputMerAligner, "%s/%s-merAlignerOutput_%d_Read%d", base_dir, libname, MYTHREAD, readNumber);
       get_rank_path(outputMerAligner, MYTHREAD);
       laneFD1 = fopen(outputMerAligner, "r");
       if (laneFD1 == NULL) {
           if (readNumber == 1) 
               DIE("Could not open %s!\n", outputMerAligner); 
           break; // otherwise okay -- may be unpaired library
       }
   
       cur_pos = 0;
#ifdef DEBUG
    int report = 0;
#endif
   
    start_io = UPC_TICKS_NOW();
   
    char *fgets_result = fgets(align, MAX_LINE_SIZE, laneFD1);
   
    end_io = UPC_TICKS_NOW();
    read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);
   
    while ( fgets_result != NULL )
    {
      
      splitRes = splitAlignment(align, cur_read, &qStart, &qStop, &qLength, &subject, &sStart, &sStop, &sLength, &strand);

      if ( splitRes == 0 ) {
         /* Found a "guard" alignment */
         /* Process the alignment iff current read alignments  >= 2 */
         if (cur_pos > 1) {
            processAlignment(&cur_alignments, cur_pos, &splints_found, outFD);
         } else if (cur_pos == 1) {
            my_sigleton++;
         }
         /* Reset the ptr in the cur_alignments data structure */
         cur_pos = 0;
      } else {
         /* If the current read (just read) is different than the previous one -> process the previous alignments */
         if (strcmp(cur_read, prev_read) != 0) {
            if (cur_pos > 1) {
               processAlignment(&cur_alignments, cur_pos, &splints_found, outFD);
            }  else if (cur_pos == 1) {
               my_sigleton++;
            }
            /* Reset the ptr in the cur_alignments data structure */
            cur_pos = 0;
            strcpy(prev_read, cur_read);
         }
         
         /* Test the status of the current read */
         testResult = findReadStatus(strand, qStart, qStop, qLength, sStart, sStop, sLength, fivePrimeWiggleRoom,threePrimeWiggleRoom, &startStatusRes, &endStatusRes, &locationRes, &my_uninformative, &my_truncated);
         
         if (testResult == SUCCESS) {
            /* Store alignment info in the data structure */
            if (cur_pos == 0) {
               strcpy(cur_alignments.read_name, cur_read);
            }
            cur_alignments.read_length = qLength;
            cur_alignments.aligned_read_locs[2*cur_pos] = qStart;
            cur_alignments.aligned_read_locs[2*cur_pos+1] = qStop;
            cur_alignments.aligned_contig_ids[cur_pos] = subject;
            cur_alignments.aligned_contig_lengths[cur_pos] = sLength;
            cur_alignments.aligned_contig_locs[2*cur_pos] = sStart;
            cur_alignments.aligned_contig_locs[2*cur_pos+1] = sStop;
            cur_alignments.aligned_strand[cur_pos] = strand;
            cur_alignments.location[cur_pos] = locationRes;
            cur_alignments.startStatus[cur_pos] = startStatusRes;
            cur_alignments.endStatus[cur_pos] = endStatusRes;
            cur_pos++;
         }
      }
      
      start_io = UPC_TICKS_NOW();
      
      fgets_result = fgets(align, MAX_LINE_SIZE, laneFD1);
      
      end_io = UPC_TICKS_NOW();
      read_io_time += UPC_TICKS_TO_SECS(end_io - start_io);

    }
   
    /* Process the last alignment left in the cur_alignment data structure (if any) */
    if (cur_pos > 1) {
       processAlignment(&cur_alignments, cur_pos, &splints_found, outFD);
    } else if (cur_pos == 1) {
       my_sigleton++;
    }
    fclose(laneFD1);
   } // iterate over files
   
   
   UPC_ATOMIC_FADD_I64(&singletons, my_sigleton);
   UPC_ATOMIC_FADD_I64(&truncated, my_truncated);
   UPC_ATOMIC_FADD_I64(&uninformatives, my_uninformative);
   UPC_ATOMIC_FADD_I64(&tot_splints, splints_found);
   
   upc_barrier;

   end = UPC_TICKS_NOW();
   if (MYTHREAD == 0) {
      printf("Splints found: %ld\n", tot_splints);
      printf("Unused alignments:\n");
      printf("UNINFORMATIVE:\t%ld\n", uninformatives);
      printf("TRUNCATED:\t%ld\n", truncated);
      printf("SINGLETON:\t%ld\n", singletons);
      printf("\nTime for computing SPLINTS : %d seconds\n" ,((int)UPC_TICKS_TO_SECS(end-start)));
      printf("\nRead I/O time is : %f seconds\n", read_io_time);
      printf("Write I/O time is : %f seconds\n", write_io_time);
      printf("Pure computation time is : %f seconds\n", (((int)UPC_TICKS_TO_SECS(end-start)))-(write_io_time + read_io_time));

      
      char countFileName[MAX_FILE_PATH];
      sprintf(countFileName,"%s-bmaMeta-splints", libname);
      get_rank_path(countFileName, -1);
      FILE *countFD = fopen_chk(countFileName, "w+" );
      fprintf(countFD, "0\t0\t0\t%d\t%s-splints\n", readLength, libname);
      fclose(countFD);
   }
   
   upc_barrier;
   free(lineBuffers);

   if (!MYTHREAD)
      printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
            ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
   return 0;
}
