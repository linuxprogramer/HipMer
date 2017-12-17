/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#include <argp.h>
#include <upc.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#include "utils.h"
#include "config.h"
#include "tracing.h"

const char *argp_program_version = "merauder 4.0 c";

#define MAX_LIBRARIES 1024
#define LIB_NAME_LEN 100

static void get_library_names(config_t *c, char *s) 
{
    c->num_lib_names = 0;
    int len = strlen(s);
    if (len >= MAX_LIBRARIES*4) 
        DIE("Too long a library string\n"); 
    int start = 0;
    int name_i = 0;
    for (int i = 0; i < len; i++) {
        if (s[i] == ',') {
            int l = i - start;
            strncpy(c->library_names[name_i], s + start, l);
            c->library_names[name_i][l] = 0;
            c->num_lib_names++;
            //serial_printf("Using lib %d: %s\n", name_i, c->library_names[name_i]);
            name_i++;
            if (name_i >= MAX_LIBRARIES) 
                DIE("Too many libraries!\n"); 
            start = i + 1;
        }
    }
    strcpy(c->library_names[name_i], s + start);
    //serial_printf("Using lib %d: %s\n", name_i, c->library_names[name_i]);
    c->num_lib_names++;
}

static int get_int_array(int max_entries, int *flags, char *s) 
{
    char *tmp = s;
    for (int i = 0; ;i++) {
        char *token = strtok(tmp, ",");
        if (!token)
            break;
        if (i >= max_entries)
            return 0;
        flags[i] = atoi(token);
        tmp = NULL;
    }
    return 1;
}


static error_t parse_opt(int key, char *optarg, struct argp_state *state) 
{
    config_t *c = state->input;
	switch (key) {
	case 'i': 
        if (!get_int_array(MAX_LIBRARIES, c->insert_size, optarg))
            DIE("Too many insert sizes\n");
        return 0;
    case 'm': c->mer_size = atoi(optarg); return 0;
    case 'D': c->min_depth = atoi(optarg); return 0;
    case 'P': c->poly_mode = 1; return 0;
    case 'v': c->verbose = atoi(optarg); return 0;
    case 'R': c->exclude_repeats = atof(optarg); return 0;
    case 'Q': c->qual_offset = atoi(optarg); return 0;
    case 'A': c->aggressive_closures = 1; return 0;
    case 'F': 
        if (!get_int_array(MAX_LIBRARIES, c->five_prime_wiggle_room, optarg))
            DIE("Too many 5p\n");
        return 0;
    case 'T': 
        if (!get_int_array(MAX_LIBRARIES, c->three_prime_wiggle_room, optarg))
            DIE("Too many 3p\n");
        return 0;
    case 'U': c->truncate = atoi(optarg); return 0;
    case 'r': 
        if (!get_int_array(MAX_LIBRARIES, c->reverse_complement, optarg))
            DIE("Too many reverse complements\n");
        return 0;
    case 'p': c->pair_projection = 0; return 0;
    case 'I': 
        if (!get_int_array(MAX_LIBRARIES, c->insert_sigma, optarg))
            DIE("Too many insert sigmas\n");
        return 0;
    case 'l':  c->max_read_len = atoi(optarg); return 0;
    case 'g': c->print_gaps = 1; return 0;
    case 'x': c->print_closures = 1; return 0;
    case 'G': c->min_gap_size = atoi(optarg); return 0;
    case 'L': c->scaff_len_cutoff = atoi(optarg); return 0;
    case 'b': get_library_names(c, optarg); return 0;
    case 's': strncpy(c->sr_file, optarg, MAX_FILE_PATH); return 0;
    case 'c': strncpy(c->contigs_file, optarg, MAX_FILE_PATH); return 0;
    case 'a': strncpy(c->mer_file, optarg, MAX_FILE_PATH); return 0;
    case 'B': strncpy(c->base_dir, optarg, MAX_FILE_PATH); return 0;
	default: return ARGP_ERR_UNKNOWN;
	}
    return ARGP_ERR_UNKNOWN;
}

static void fail_reqd(struct argp *ap, const char *msg) 
{
    argp_help(ap, stdout, ARGP_HELP_USAGE, (char *)msg);
    upc_global_exit(1);
}

static void set_argp_option(struct argp_option *opt, int key, const char *arg, const char *doc) 
{
    opt->name = 0;
    opt->key = key;
    opt->arg = arg;
    opt->flags = 0;
    opt->doc = doc;
    opt->group = 0;
}

static void tprint_array(char *label, int *array, int num) 
{
    tprintf("%s", label);
    for (int i = 0; i < num; i++) 
        tprintf("%d%c", array[i], i == num - 1 ? ' ' : ',');
    tprintf("\n");
}

static shared int gaps_file_exists;

void get_config(int argc, char** argv, config_t *cfg) 
{
    shared config_t *allcfg = upc_all_alloc(1, sizeof(config_t));
    upc_barrier;

    if (MYTHREAD == 0) {
        for (int i = 0; i < argc; i++)
            serial_printf("%s ", argv[i]);
        serial_printf("\n");
        struct argp_option options[40];
        int i = 0;
        // need this ugliness because of UPC
        set_argp_option(&options[i++], 'i', "char*", "Insert size per library (comma separated)");
        set_argp_option(&options[i++], 'I', "char*", "Insert sigma per library (comma separated)");
        set_argp_option(&options[i++], 'm', "int", "Mer size");
        set_argp_option(&options[i++], 'D', "int", "Minimum depth (default 2)");
        set_argp_option(&options[i++], 'P', 0, "Polymorphic mode");
        set_argp_option(&options[i++], 'v', "int", "Verbose mode, 1, 2 or 3");
        set_argp_option(&options[i++], 'R', "double", "Exclude repeats");
        set_argp_option(&options[i++], 'Q', "int", "Quality offset (default 33)");
        set_argp_option(&options[i++], 'A', 0, "Aggressive closures");
        set_argp_option(&options[i++], 'd', "char*", "Data directory");
        set_argp_option(&options[i++], 'f', "char*", "FASTQ sequence file glob");
        set_argp_option(&options[i++], 's', "char*", "Scaffold report file pattern");
        set_argp_option(&options[i++], 'c', "char*", "Contigs FASTA file pattern");
        set_argp_option(&options[i++], 'a', "char*", "Meraligner file pattern");
        set_argp_option(&options[i++], 'b', "char*", "Meraligner file indexes (comma separated)");
        set_argp_option(&options[i++], 'B', "char*", "Base directory for temporary files");
        set_argp_option(&options[i++], 'F', "char*", "5' wiggle room per library (comma separated)");
        set_argp_option(&options[i++], 'T', "char*", "3' wiggle room per library (comma separated)");
        set_argp_option(&options[i++], 'l', "int", "Max read length over all libraries");
        set_argp_option(&options[i++], 'U', "int", "Truncate");
        set_argp_option(&options[i++], 'r', "char*", "Reverse complement per library (comma separated) ");
        set_argp_option(&options[i++], 'p', 0, "No pair projection");
        set_argp_option(&options[i++], 'G', "int", "Min gap size");
        set_argp_option(&options[i++], 'L', "int", "Scaffold length cutoff");
        set_argp_option(&options[i++], 'g', 0, "Print gaps");
        set_argp_option(&options[i++], 'x', 0, "Print closures");
        set_argp_option(&options[i], 0, 0, 0);

        config_t mycfg;
        // defaults
        mycfg.min_depth = 2;
        mycfg.poly_mode = 0;
        mycfg.verbose = 0;
        mycfg.exclude_repeats = 0;
        mycfg.qual_offset = 33;
        mycfg.aggressive_closures = 0;
        for (int j = 0; j < MAX_LIBRARIES; j++) {
            mycfg.five_prime_wiggle_room[j] = 5;
            mycfg.three_prime_wiggle_room[j] = 5;
        }
        mycfg.max_read_len = 0;
        mycfg.truncate = 0;
        memset(mycfg.reverse_complement, 0, sizeof(mycfg.reverse_complement));
        mycfg.pair_projection = 1;
        mycfg.print_gaps = 0;
        mycfg.print_closures = 0;
        mycfg.min_gap_size = 10;
        mycfg.scaff_len_cutoff = 1000;
        mycfg.num_lib_names = 0;
        // required 
        memset(mycfg.insert_size, 0, sizeof(mycfg.insert_size));
        memset(mycfg.insert_sigma, 0, sizeof(mycfg.insert_sigma));
        mycfg.mer_size = -1;
		strcpy(mycfg.sr_file, "");
        strcpy(mycfg.contigs_file, "");
        strcpy(mycfg.mer_file, "");
        strcpy(mycfg.base_dir, ".");

        static struct argp ap;
		ap.options = options;
		ap.parser = parse_opt;
		ap.args_doc = "";
		ap.doc = "Options:";

        int arg_index = 1;
        int argParsed = argp_parse(&ap, argc, argv, ARGP_IN_ORDER, &arg_index, &mycfg);
        if (argParsed != 0)  {
            fprintf(stderr, "argp_parse returned %d\n", argParsed);
            fail_reqd(&ap, "Failed to parse the command line options\n");
        }
        mycfg.overall_max_insert_size = 0;
        for (int l = 0; l < mycfg.num_lib_names; l++) {
            if (!mycfg.insert_size[l])
                fail_reqd(&ap, "insert size is required for each lib (-i int,int,...)\n");
            if (!mycfg.insert_sigma[l] && !gaps_file_exists)
                fail_reqd(&ap, "insert sigma is required for each lib (-I int,int,...)\n");
            if (mycfg.insert_size[l] > mycfg.overall_max_insert_size)
                mycfg.overall_max_insert_size = mycfg.insert_size[l];
            // if the wiggles were set to 0, that means they are not set and should be at the default of 5
            if (!mycfg.five_prime_wiggle_room[l])
                mycfg.five_prime_wiggle_room[l] = 5;
            if (!mycfg.three_prime_wiggle_room[l])
                mycfg.three_prime_wiggle_room[l] = 5;
        }
        // default max read lens to 101
        if (!mycfg.max_read_len)
            fail_reqd(&ap, "max read len is required (-l int)\n");
        if (mycfg.mer_size == -1) 
            fail_reqd(&ap, "mer size is required (-m int)\n");
		if (mycfg.sr_file[0] == 0)
			fail_reqd(&ap, "need to specify scaffold report file (-s char*)\n");
		if (mycfg.contigs_file[0] == 0)
			fail_reqd(&ap, "need to specify contigs FASTA file (-c char*)\n");
		if (mycfg.mer_file[0] == 0)
			fail_reqd(&ap, "need to specify Meraligner file (-a char*)\n");
        
        upc_memput(allcfg, &mycfg, sizeof(config_t));
    }
    upc_barrier;
    upc_memget(cfg, allcfg, sizeof(config_t));
    set_output_dir(cfg->base_dir);
    tprintf("%s\n", argp_program_version);
    tprint_array("  Insert_size:          ", cfg->insert_size, cfg->num_lib_names);
    tprint_array("  Insert sigma:         ", cfg->insert_sigma, cfg->num_lib_names);
    tprintf("  mer_size:             %d\n", cfg->mer_size);
    tprintf("  min_depth:            %d\n", cfg->min_depth);
    tprintf("  qual_offset:          %d\n", cfg->qual_offset);
    tprintf("  exclude repeats:      %.4f\n", cfg->exclude_repeats);
    tprintf("  polymorphic mode:     %s\n", cfg->poly_mode ? "true" : "false");
    tprintf("  verbose mode:         %d\n", cfg->verbose);
    tprintf("  aggressive closures:  %s\n", cfg->aggressive_closures ? "true" : "false");
    tprint_array("  5' wiggle room:       ", cfg->five_prime_wiggle_room, cfg->num_lib_names);
    tprint_array("  3' wiggle room:       ", cfg->three_prime_wiggle_room, cfg->num_lib_names);
    tprintf("  Max read len:         %d\n", cfg->max_read_len);
    tprintf("  Truncate:             %d\n", cfg->truncate);
    tprint_array("  Reverse complement:   ", cfg->reverse_complement, cfg->num_lib_names);
    tprintf("  Pair projection:      %s\n", cfg->pair_projection ? "true" : "false");
    tprintf("  Print gaps:           %s\n", cfg->print_gaps ? "true" : "false");
    tprintf("  Print closures:       %s\n", cfg->print_closures ? "true" : "false");
    tprintf("  Min gap size:         %d\n", cfg->min_gap_size);
    tprintf("  Scaffold len cutoff:  %d\n", cfg->scaff_len_cutoff);
    tprintf("  Scaffold report file: %s\n", cfg->sr_file);
	tprintf("  Contigs FASTA file:   %s\n", cfg->contigs_file);
    tprintf("  Meraligner file:      %s\n", cfg->mer_file);
    tprintf("  Base dir:             %s\n", cfg->base_dir);
	tprintf("  library names:        ");
	for (int i = 0; i < cfg->num_lib_names; i++) 
        tprintf("%s%c", cfg->library_names[i], i == cfg->num_lib_names - 1 ? ' ' : ',');
	tprintf("\n");
    tflush();
}

