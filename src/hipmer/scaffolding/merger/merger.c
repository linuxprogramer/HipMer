#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "../../common/optlist.h"
#include "../../common/common.h"
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

#define MAX_LINE_SIZE 1000
#define MAX_FILES 10

int main(int argc, char **argv) {
   upc_tick_t start_time = upc_ticks_now();

   int nFiles;
   char line[MAX_LINE_SIZE+1];
   char temp_name[MAX_FILE_PATH+1];
   
   option_t *optList, *thisOpt;
   optList = NULL;
   optList = GetOptList(argc, argv, "P:f:s:");
   char *fileName, *suffix = NULL;
   
   /* Process the input arguments */
   while (optList != NULL) {
      thisOpt = optList;
      optList = optList->next;
      switch (thisOpt->option) {
         case 'P':
            nFiles = atoi(thisOpt->argument);
            break;
         case 'f':
            fileName = thisOpt->argument;
            break;
         case 's':
            suffix = thisOpt->argument;
            break;
         default:
            break;
      }
      free(thisOpt);
   }

   serial_printf("Merging %s%s files from %d to %d\n", fileName, suffix == NULL ? "" : suffix, nFiles, THREADS);

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
   int64_t *fileSizes = (int64_t*) malloc_chk(myFiles*sizeof(int64_t));
   FILE **FD_array = (FILE **) malloc_chk(myFiles*sizeof(FILE*));
   
   int fd;

   for (i = my_offset; i < my_offset + myFiles; i++) {
       /* Adapt the name to the appropriate THREAD id */
       sprintf(temp_name, "%s_%d%s", fileName, i, suffix == NULL ? "" : suffix);
       curReadFD = fopen_rank_path(temp_name, "r", i);
       fd = fileno(curReadFD);
       if(fstat(fd,&fileStat) < 0) return 1;
       cur_file_size = fileStat.st_size;
       fileSizes[k] = cur_file_size;
       FD_array[k] = curReadFD;
       k++;
   }
   
   int64_t maxFileSize = 0;
   for (j = 0; j < myFiles; j++) {
      if (fileSizes[j] > maxFileSize) {
         maxFileSize = fileSizes[j];
      }
   }
   
   /* Alocate buffer of maximum size */
   char *inputBuffer = (char*) malloc((maxFileSize+1) * sizeof(char));

   k = 0;
   
   sprintf(temp_name, "merged_%s_%d", fileName, MYTHREAD);
   curWriteFD = fopen_rank_path(temp_name, "w", MYTHREAD);
   for (i = my_offset; i < my_offset + myFiles; i++) {
       /* Adapt the name to the appropriate THREAD id */
       curReadFD = FD_array[k];
       if (fread(inputBuffer, 1, fileSizes[k], curReadFD))
           fwrite(inputBuffer, 1, fileSizes[k], curWriteFD);
       fclose(curReadFD);
       k++;
   }
   fclose(curWriteFD);
   
   upc_barrier;
   
   if (!MYTHREAD)
      printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
            ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));

   return 0;
}
