#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <string.h>
#include <upc.h>
#include <upc_tick.h>
#include <libgen.h>

#include "../fqreader/fq_reader.h"
#include "../common/optlist.h"

#define MAX_SEQ_FILES 12

#define MAX_LIBRARIES 1024
#define LIB_NAME_LEN 100
#define MAX_FILE_PATH 384 

static char library_names[MAX_LIBRARIES][LIB_NAME_LEN];
static int num_lib_names = 0;

static void get_library_names(char *s) 
{
    num_lib_names = 0;
    int len = strlen(s);
    if (len >= MAX_LIBRARIES*4) 
        DIE("Too long a library string\n"); 
    int start = 0;
    int name_i = 0;
    for (int i = 0; i < len; i++) {
        if (s[i] == ',') {
            int l = i - start;
            strncpy(library_names[name_i], s + start, l);
            library_names[name_i][l] = 0;
            num_lib_names++;
            serial_printf("Using lib %d: %s\n", name_i, library_names[name_i]);
            name_i++;
            if (name_i >= MAX_LIBRARIES) 
                DIE("Too many libraries!\n"); 
            start = i + 1;
        }
    }
    strcpy(library_names[name_i], s + start);
    serial_printf("Using lib %d: %s\n", name_i, library_names[name_i]);
    num_lib_names++;
}

int main(int argc, char** argv)
{
    option_t *this_opt;
    option_t *opt_list = GetOptList(argc, argv, "f:");
    char *libnames_ext = NULL;
    while (opt_list) {
        this_opt = opt_list;
        opt_list = opt_list->next;
        switch (this_opt->option) {
        case 'f':
            libnames_ext = this_opt->argument;
            break;
        default:
            serial_printf("Invalid option %c\n", this_opt->option);
            upc_global_exit(1);
      }
    }
    if (!libnames_ext) {
        serial_printf("Usage: %s -f fofn_base_name\n", argv[0]);
        upc_global_exit(1);
    }
    serial_printf("Loading FASTQ files from %s\n", libnames_ext);
    get_library_names(libnames_ext);

    upc_tick_t t = upc_ticks_now();

    for (int libi = 0; libi < num_lib_names; libi++) {
		char fofn[MAX_FILE_PATH], fname[MAX_FILE_PATH];
		snprintf(fofn, MAX_FILE_PATH, "%s.fofn", library_names[libi]);
		serial_printf("Processing fastq files from %s...\n", fofn);
		FILE *fofn_fd = fopen_chk(fofn, "r");
		while (fgets(fname, MAX_FILE_PATH, fofn_fd)) { 
			fname[strlen(fname)-1] = '\0';
            serial_printf("Loading %s into memory on %d threads\n", fname, THREADS);
            fq_reader_t fqr = create_fq_reader();
            int err = load_fq(fqr, fname);
            if (err < 0)
                printf("Could not load %s into memory: %s\n", fname, strerror(-err));
        }
    }
	upc_barrier;
	serial_printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
				  ((double)upc_ticks_to_ns(upc_ticks_now() - t) / 1000000000.0));
    return 0;
}
