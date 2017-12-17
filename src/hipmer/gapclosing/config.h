/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#ifndef __CONFIG_H
#define __CONFIG_H

#include <limits.h>
#include "../common/common.h"
#include "utils.h"

#define MAX_LIBRARIES 1024
#define LIB_NAME_LEN 100

typedef struct {
    int insert_size[MAX_LIBRARIES];
    int insert_sigma[MAX_LIBRARIES];
    int overall_max_insert_size;
    int mer_size;
    int min_depth;
    int poly_mode;
    int verbose;
    double exclude_repeats;
    int qual_offset;
    int aggressive_closures;
    int five_prime_wiggle_room[MAX_LIBRARIES];
    int three_prime_wiggle_room[MAX_LIBRARIES];
    int max_read_len;
    int truncate;
    int reverse_complement[MAX_LIBRARIES];
    int pair_projection;
    int print_gaps;
    int print_closures;
    int min_gap_size;
    int scaff_len_cutoff;
    char library_names[MAX_LIBRARIES][LIB_NAME_LEN];
    int num_lib_names;
    char sr_file[MAX_FILE_PATH + 1];
    char contigs_file[MAX_FILE_PATH + 1];
    char mer_file[MAX_FILE_PATH + 1];
    char base_dir[MAX_FILE_PATH + 1];
} config_t;

void get_config(int argc, char** argv, config_t *cfg);

#endif
