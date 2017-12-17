#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <upc.h>
#include <upc_tick.h>
#include <string.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <libgen.h>
#include "../../common/optlist.h"
#include "../../common/common.h"

#define MAX_LINE_SIZE 1000
#define MAX_FILES 10

int main(int argc, char **argv) {
   upc_tick_t start_time = upc_ticks_now();

   int nFiles;
   char line[MAX_LINE_SIZE+1];
   char cpy_line[MAX_LINE_SIZE+1];
   char temp_name[MAX_LINE_SIZE+1];
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "P:f:s:");
   char *fileNames;
   char *suffixName;
   
   /* Process the input arguments */
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'P':
            nFiles = atoi(thisOpt->argument);
            break;
         case 's':
            suffixName = thisOpt->argument;
            break;
         case 'f':
            fileNames = thisOpt->argument;
            break;
         default:
            break;
      }
      free(thisOpt);
   }
   
   serial_printf("Splitting %s files from %d to %d\n", fileNames, THREADS, nFiles);

   int cur_procs = THREADS;
   int remainder = nFiles % cur_procs;
   int teamFiles = (nFiles - remainder) / cur_procs;
   int *filesPerTeam = (int*) malloc(THREADS * sizeof(int));
   int i;
   
   for (i=0; i<THREADS; i++) {
      filesPerTeam[i] = teamFiles;
      if (i < remainder) {
         filesPerTeam[i]++;
      }
   }
   
   int myFiles = filesPerTeam[MYTHREAD];
   
   int my_offset = 0;
   for (i=0; i<MYTHREAD; i++) {
      my_offset += filesPerTeam[i];
   }
   
   struct stat fileStat;
   int j;
   int64_t cur_buf_size = 0;
   int64_t cur_file_size = 0;
   FILE *curReadFD, *curWriteFD;
   
   int k = 0;
   FILE **FD_array = (FILE **) malloc(myFiles*sizeof(FILE*));
   
   for (i = my_offset; i < my_offset + myFiles; i++) {
      /* Adapt the name to the appropriate THREAD id */
      sprintf(temp_name, "%s_%s_%d", fileNames, suffixName, i);
      curWriteFD = fopen_rank_path(temp_name, "a+", i);
      FD_array[k] = curWriteFD;
      k++;
   }

   sprintf(temp_name, "%s_%d", fileNames, MYTHREAD);
   curReadFD = fopen_rank_path(temp_name, "r", MYTHREAD);
   
   k = -1;
   char *token;
   char *aux;
   int scaffoldId = 0;
   int prevScaffoldId = -1;
   
   char *fgets_result = fgets(line, MAX_LINE_SIZE, curReadFD);
   while ( fgets_result  != NULL) {
      strcpy(cpy_line, line);
      token = strtok_r(line, "\t", &aux);
      scaffoldId = atoi(token+8);
      
      if (scaffoldId != prevScaffoldId) {
         k++;
         k = k % myFiles;
         curWriteFD = FD_array[k];
         prevScaffoldId = scaffoldId;
      }
      
      fprintf(curWriteFD, "%s", cpy_line);
      fgets_result = fgets(line, MAX_LINE_SIZE, curReadFD);
   }
   
   fclose(curReadFD);
   for (k = 0; k<myFiles; k++) {
      fclose(FD_array[k]);
   }

   upc_barrier;
   
   if (!MYTHREAD)
      printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
            ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));

   return 0;
}
