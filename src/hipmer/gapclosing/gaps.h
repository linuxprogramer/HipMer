/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#ifndef __GAPS_H
#define __GAPS_H

#include <stdio.h>
#include "htable.h"
#include "config.h"

#include "../common/Buffer.h"

typedef struct {
    int id;
    double length;
    double depth;
    double weight;
    int ncontigs;
    shared [] int *contig_ids;
} scaffold_t;

typedef struct {
    int id;
    shared [] char *seq;
    int seq_len;
    char strand;
    int gapi_3;
    int gapi_5;
    int s_start;
    int s_end;
} contig_t;

typedef struct {
    int id;
    int scaff_id;
    int contig1_id;
    int contig2_id;
    int contig_ext1;
    int contig_ext2;
    char *primer1;
    char *primer2;
    int size;
    double uncertainty;
    int start;
    int end;
    int max_reads;
    int nreads;
    int len_closure;
    shared [] int *reads_lens;
    shared [] char *reads_nts;
    shared [] char *reads_quals;
    shared [] char *closure;
} gap_t;

typedef struct {
    shared [] gap_t *gaps;
    int ngaps;
    int gap_start;
    int gap_end;
} gap_set_t;

#define NULL_ON_FAIL 0
#define ABORT_ON_FAIL 1

void init_globals(config_t *cfg, shared gap_set_t *gaps);
shared [] gap_t *get_gap(int gi, int line, const char *fname, int abort_on_fail);
char *get_gap_info(gap_t *gap, char *buf, int nchar);
char *get_report_line(gap_t *gap, Buffer buf);
void fprint_gap(FILE *f, gap_t *gap);
void process_contigs_file(void);
void process_scaffold_report_file(void);
void process_meraligner_file(char *library_name, int reverse_complement, int insert_size, 
                             int insert_sigma, int five_prime_wiggle, int three_prime_wiggle);
void process_fastq_file(char *library_name);
void print_gaps(void);
shared contig_t *get_contig(int id);
#define get_contig_chk(contig, id)                      \
    do {                                                \
        contig = get_contig(id);                        \
        if (!contig)                                    \
            DIE("Could not get contig for %d\n", id);   \
    } while (0);
shared scaffold_t *get_scaffold(int id);
char *get_primer(shared contig_t *contig, int contig_ext, int reverse_flag);
shared contig_t *get_contigs(void);
shared scaffold_t *get_scaffolds(void);
int get_max_contigs(void);
int get_max_scaffolds(void);
void load_balance_gaps(void);

#endif
