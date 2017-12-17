/*
  Calculates the N50 value from a collection of scaffold files (srf).
  Each srf file has contig entries of the form:

  Scaffold<si>  CONTIG<sc> <+|->Contig<numbers>  start_contig  end_contig  depth

  All we need is the end_contig of the last contig for each scaffold. 
  Then we use those to compute the N50. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <upc.h>
#include <upc_tick.h>

#include "../../common/optlist.h"
#include "../../common/common.h"

#define MAX_LINE_SIZE 255


static int cmp_int(const void *a, const void *b)
{
    return (*(int*)a - *(int*)b);
}

int main(int argc, char **argv) 
{
   upc_tick_t start_time = upc_ticks_now();
   int num_scaffolds = 0;
   char *file_prefix = NULL;
   option_t *opt_list = GetOptList(argc, argv, "n:f:");
   while (opt_list != NULL) {
       option_t *this_opt = opt_list;
       opt_list = opt_list->next;
       switch (this_opt->option) {
       case 'n':
           num_scaffolds = atoi(this_opt->argument);
           break;
       case 'f':
           file_prefix = this_opt->argument;
           break;
         default:
            break;
      }
      free(this_opt);
   }
   serial_printf("Running on files %s, with %d scaffolds\n", file_prefix, num_scaffolds);
   if (!num_scaffolds || !file_prefix) {
       serial_printf("Usage: %s -n num_scaffolds -f file_prefix\n", argv[0]);
       upc_global_exit(0);
   }
   
   shared int *scaff_lens = upc_all_alloc(num_scaffolds, sizeof(int));

   int scaff_word_len = strlen("Scaffold");
   char fname[MAX_FILE_PATH];
   sprintf(fname, "%s_%d", file_prefix, MYTHREAD);
   FILE *f = fopen_chk(fname, "r");
   char line[MAX_LINE_SIZE];
   int scaff_id = -1;
   int prev_scaff_id = -1;
   int num_scaffs = 0;
   int scaff_len = 0;
   int prev_scaff_len = 0;
   while (fgets(line, MAX_LINE_SIZE, f)) {
       char *aux;
       char *token = strtok_r(line, "\t", &aux);
       if (strncmp(token, "Scaffold", scaff_word_len) != 0)
           DIE("Invalid format in SRF file, should start with 'Scaffold<N>': %s\n", line);
       prev_scaff_id = scaff_id;
       scaff_id = atoi(token + scaff_word_len);
       if (scaff_id >= num_scaffolds)
           DIE("Too many scaffolds for -n parameter %d\n", num_scaffolds);
       token = strtok_r(NULL, "\t", &aux);
       if (token[0] != 'C')
           continue;
       token = strtok_r(NULL, "\t", &aux);       
       token = strtok_r(NULL, "\t", &aux);
       token = strtok_r(NULL, "\t", &aux);
       prev_scaff_len = scaff_len;
       scaff_len = atoi(token);
       if (scaff_id != prev_scaff_id && prev_scaff_id != -1) {
           // note that the scaff id should never be repeated across files
           if (scaff_lens[prev_scaff_id])
               DIE("Duplicate scaffold id found %d\n", prev_scaff_id);
           scaff_lens[prev_scaff_id] = prev_scaff_len;
       }
   }
   // get the last one
   scaff_lens[scaff_id] = scaff_len;
   fclose(f);
   upc_barrier;
   if (!MYTHREAD) {
       int *local_scaff_lens = malloc_chk(num_scaffolds * sizeof(int));
       int tot_len = 0;
       for (int i = 0; i < num_scaffolds; i++) {
           local_scaff_lens[i] = scaff_lens[i];
           tot_len += local_scaff_lens[i];
       }
       qsort(local_scaff_lens, num_scaffolds, sizeof(int), cmp_int);
       //printf("Scaff lens:\n");
       //for (int i = 0; i < num_scaffolds; i++) {
       //    printf("%d %d\n", i, local_scaff_lens[i]);
       //}
       // compute the N50
       long running_tot = 0;
       double idx = 0.5;
       for (int i = num_scaffolds - 1; i >= 0; i--) {
           running_tot += local_scaff_lens[i];
           if (running_tot >= idx * tot_len) {
               int nx = (int)(idx * 100);
               serial_printf("\tN%d %d L%d %d\n", nx, local_scaff_lens[i], nx, i + 1);
               //printf("%d\n", local_scaff_lens[i]);
               break;
           }
       }
   }
   serial_printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
                 ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
   return 0;
}


