/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 *
 * Based on merauder.pl by Jarrod Chapman <jchapman@lbl.gov> Tue Jun 2
 * 07:48:59 PDT 2009 Copyright 2009 Jarrod Chapman. All rights reserved.
*/

/*
 * Three files are read in:
 * 1. The gap data, where all of the info is stored in an array of gap structs. 
 * 2. The contigs from a FASTA file, where only the contigs for relevant
 * sequences are stored.
 * 3. The scaffold report file, which is used to determine the start, end and
 * depth of the scaffold.
 *
 * Each thread outputs to a separate file.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <upc.h>
#include <upc_collective.h>
#include <ctype.h>
#include <upc_tick.h>
#include <libgen.h>
#include <unistd.h>

#include "htable.h"
#include "gaps.h"
#include "utils.h"
#include "config.h"
#include "timers.h"
#include "tracing.h"
#include "../common/Buffer.h"

// for debugging
//#define SOLO_GAP 62

//#define DBG_SCAFF_FASTA

#define CHECK_SEQ_POS(x) do {                               \
        char label[500];                                    \
        sprintf(label, "  in %s:%d\n", __FILE__, __LINE__); \
        check_seqs(x, label);                               \
    } while(0);

#define MAX_SCAFF_SEQ_LEN 1000000
#define FASTA_LINE_WIDTH 50

#define tprintf_vvv(fmt,...) if (_cfg.verbose > 2) tprintf(fmt, ##__VA_ARGS__)
#define tprintf_vv(fmt,...) if (_cfg.verbose >= 2) tprintf(fmt, ##__VA_ARGS__)
#define tprintf_v(fmt,...) if (_cfg.verbose == 1) tprintf(fmt, ##__VA_ARGS__)

#define DIRN_RIGHT 0
#define DIRN_LEFT 1

static const char *DIRNS[2] = {"right", "left"};

// Limit for all htables. Exceeding this limit won't cause a failure, it will
// just cause inefficiencies, since collisions are resolved through linked lists
// note that this number gets rounded up to a prime close to 2^i
//#define MAX_HTABLE_ENTRIES 60000
#define MAX_HTABLE_ENTRIES 100000
// For output only
#define MAX_BUF 2000
#define MAX_WALK_LEN 8000
#define MAX_BEST_GUESS_SEQ_LEN 2000

#define GET_READ(s, i) ((s) + ((_cfg.max_read_len + 1) * (i)))

typedef struct {
    char *span;
    int freq;
} span_info_t;

typedef struct {
    char *mer;
    int freqs[4][4];
} qual_freqs_t;

typedef struct {
    char base;
    int rating;
    int n_hi_q;
    int n_ext;
} qual_t;

typedef struct {
    char *mer;
    char base;
    char fork_bases[4];
    int fork_n_ext[4];
} mer_base_t;

typedef struct {
    shared int *scaff_lens;
    shared int *contig_lens;
    int scaff_tot;
    long scaff_len_tot;
    int contig_tot;
    long contig_len_tot;
} assembly_stats_t;

static config_t _cfg;

static int _nwalks = 0;
static int _max_htable_entries = MAX_HTABLE_ENTRIES;

static int _nsuccess = 0;
static int _nfailure = 0;

DEFINE_HTABLE_TYPE(int32)
DEFINE_HTABLE_TYPE(span_info);
DEFINE_HTABLE_TYPE(qual_freqs);
DEFINE_HTABLE_TYPE(mer_base);

static const char BASES[4] = {'A', 'C', 'G', 'T'};

static inline char *get_seq(shared contig_t *contig, char *primer, int dirn)
{
    if (!contig->seq)
        return NULL;
    char *seq = malloc_chk(contig->seq_len + 1);
    UPC_MEMGET_STR(seq, contig->seq, contig->seq_len + 1);
    if (contig->strand == '-') 
        switch_code(reverse(seq));
#ifdef CONFIG_SANITY_CHECK
    if ((dirn == DIRN_LEFT && endswith(seq, primer) != 0) ||
        (dirn == DIRN_RIGHT && strncmp(seq, primer, strlen(primer)) != 0)) {
        DIE("Contig sequence doesn't match primer at %s (strand %c):\nseq (%d): '%s'\nprimer: '%s'\n", 
             dirn == DIRN_LEFT ? "end" : "start", contig->strand, contig->seq_len, seq, primer);
    }
#endif
    return seq;
}

static char *last_strstr(char *haystack, char *needle)
{
    char *temp = haystack, *before = 0;
    while ((temp = strstr(temp, needle))) 
        before = temp++;
    return before;
}

static char *span(char *primer1, char *primer2, int *reads_lens, char *reads_nts, char *reads_quals,
                  int nreads) 
{
	tprintf_vvv("Attempting span: %s -> %s (%d read(s) available)\n", primer1, primer2, nreads);

    char *ret = NULL;
    htable_t pure_spans = create_htable(MAX_HTABLE_ENTRIES, "pure_spans");
    span_info_t *first_span = NULL;
    char *span = NULL;
    char *nts, *p1_substr, *p2_substr;
    span_info_t *span_info;
    for (int i = 0; i < nreads; i++) {
        nts = GET_READ(reads_nts, i);
        if (!nts)
            DIE("not nts, length %d\n", reads_lens[i]);
        p1_substr = strstr(nts, primer1);
        // get the longest span
        p2_substr = last_strstr(nts, primer2);
        if (p1_substr && p2_substr) {
            if (p1_substr <= p2_substr) {
                // a string starting with primer1 and finishing with primer2
                span = strdup(p1_substr);
                span[p2_substr + strlen(primer2) - p1_substr] = '\0';
                span_info = htable_get_span_info(pure_spans, span);
                if (span_info) {
                    span_info->freq++;
                    free(span);
                } else {
                    span_info = malloc_chk(sizeof(span_info_t));
                    span_info->freq = 1;
                    span_info->span = span; 
                    CHECK_ERR(htable_put_span_info(pure_spans, span_info->span, span_info, NO_CHECK_DUPS));
                    if (!first_span)
                        first_span = span_info;
                }
            }
        }
    }
    int nspans = htable_num_entries(pure_spans);
	tprintf_vvv("%d distinct spanning sequence(s) found\n", nspans);
    if (nspans == 1 && first_span->freq > 1) {
        tprintf_vv("Unique spanning sequence found: %d/%d reads span the gap.\n",
				   first_span->freq, nreads);
        ret = strdup(first_span->span);
    } else if (_cfg.poly_mode) {
        int max_span_freq = 0;
        htable_iter_t iter = htable_get_iter(pure_spans);
        span_info_t *span_info = NULL;
        int new_len, prev_len;
        while ((span_info = htable_get_next_span_info(pure_spans, iter)) != NULL) {
            if (span_info->freq > max_span_freq) {
                span = span_info->span;
                max_span_freq = span_info->freq;
            } else if (span_info->freq == max_span_freq) {
                new_len = strlen(span_info->span);
                prev_len = strlen(span);
                if (new_len > prev_len) {
                    // always prefer the longest span for the same frequency
                    span = span_info->span;
                } else if (new_len == prev_len && strcmp(span_info->span, span) < 0) {
                    // prefer lexicographically smallest f or consistency
                    span = span_info->span;
                }
            }
        }
        free(iter);
        if (max_span_freq > 1) {
            tprintf_vv("Maximum frequency spanning sequence found: %d/%d reads span the gap.\n",
					   max_span_freq, nreads);
            ret = span;
        } else {
			tprintf_vvv("No spanning sequence with frequency > 1 found\n");
        }
    } else {
		tprintf_vvv("No unique spanning sequence with frequency > 1 found\n");
    }
    destroy_htable(pure_spans, FREE_VALS);
    return ret;
}

static int check_closure(char *closure, char *primer1, char *primer2, int gap_size, 
                         double gap_uncertainty)
{
    int bad_closure = 0;
    int p1_match = strncmp(closure, primer1, strlen(primer1));
    int p2_match = endswith(closure, primer2);

    if (p1_match != 0 || p2_match != 0) {
        bad_closure++;
		tprintf_vvv("closure [%s] rejected because it disagrees with primers\n", closure);
    }
    int closed_gap_size = strlen(closure) - 2 * _cfg.mer_size;
    int gap_diff = closed_gap_size - gap_size;

    if (abs(gap_diff) > gap_uncertainty) {
        bad_closure += 2;
		tprintf_vvv("closure [%s] rejected due to gap estimate differential: |%d| > %.4f\n",
                    closure, gap_diff, gap_uncertainty);
    }
    // In aggressive mode ignore size differences between closure and estimate
    if (_cfg.aggressive_closures && (bad_closure == 2) && (gap_size < 2 * _cfg.overall_max_insert_size)) {
        bad_closure = 0;
		tprintf_vvv("closure [%s] allowed via AGRESSIVE mode\n", closure);
    }
    return bad_closure;
}

#define add_to_report(report_note, fmt, ...) printfBuffer(report_note, fmt, __VA_ARGS__)
static void old_add_to_report(char *report_note, const char *fmt, ...) 
{
    va_list args;
    va_start(args, fmt);
    char buf[MAX_BUF];
    vsnprintf(buf, MAX_BUF - 1, fmt, args);
    va_end(args);
    strcat(report_note, buf);
}

static char *splinting_reads(gap_t *gap, int *span_check, Buffer report_note) 
{   
    START_TIMER(T_SPLINTING);
    char *span_closure = span(gap->primer1, gap->primer2, (int*)gap->reads_lens, 
                              (char*)gap->reads_nts, (char*)gap->reads_quals, gap->nreads);
    *span_check = 0;
    if (span_closure) {
        *span_check = check_closure(span_closure, gap->primer1, gap->primer2, gap->size,
                                    gap->uncertainty);
        add_to_report(report_note, "spanClosure=%s;spanCheck=%d;", span_closure, *span_check);
        if (*span_check != 0) 
            span_closure = NULL;
    }
    stop_timer(T_SPLINTING);
    return span_closure;
}

static inline int get_base_index(char base) 
{
    switch (base) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    }
    return -1;
}

static int get_valid_base_str(char *s, int max_len, int *start, int *len)
{
    *start = 0;
    int first_N = 0;
    for (int j = 0; j < *len; j++) {
        if (s[j] == 'N') {
            first_N = j;
            if (j - *start >= max_len) {
                *len = j - *start;
                break;
            } else {
                while (s[j] == 'N')
                    j++;
                if (j < *len - 1)
                    *start = j;
                else
                    *len = first_N - *start;
            }
        }
        if (j == *len - 1)
            *len -= *start;
    }
    if (*len >= max_len)
        return 1;
    else 
        return 0;
}

static int categorize_extension(int *ext)
{
    int min_viable = 3;
    if (min_viable > _cfg.min_depth) {
        tprintf_vv("Warning: in categorizeExtension minViable reset to match minDepth (%d)\n", 
				   _cfg.min_depth);
        min_viable = _cfg.min_depth;
    }

    // 0 = No votes
    // 1 = One vote
    // 2 = nVotes < minViable
    // 3 = minDepth > nVotes >= minViable, nHiQ < minViable
    // 4 = minDepth > nVotes >= minViable ; nHiQ >= minViable
    // 5 = nVotes >= minDepth ; nHiQ < minViable
    // 6 = nVotes >= minDepth ; minViable < nHiQ < minDepth 
    // 7 = nHiQ >= minDepth 

    //    Ignore q<10 bases
    //    my $n = $ext[0]+$ext[1]+$ext[2]+$ext[3];
    int n = ext[1] + ext[2] + ext[3];
    int n_hi_q = ext[2] + ext[3];
    int category = -1;
    if (n == 0) {
        category = 0;
    } else if (n == 1) {
        category = 1;
    } else if (n < min_viable) {
        category = 2;
    } else {
        if ((n < _cfg.min_depth) || (n == min_viable)) {
            if (n_hi_q < min_viable) 
                category = 3;
            else 
                category = 4;
        } else {
            if (n_hi_q < min_viable) 
                category = 5;
            else if (n_hi_q < _cfg.min_depth) 
                category = 6;
            else 
                category = 7;
        }
    }

    if (category == -1)
        DIE("Undefined extension category");
    return category;
}

static inline int cmp_quals(qual_t *q1, qual_t *q2)
{
    if (q1->rating > q2->rating)
        return 1;
    if (q1->rating < q2->rating)
        return -1;
    if (q1->n_hi_q > q2->n_hi_q)
        return 1;
    if (q1->n_hi_q < q2->n_hi_q)
        return -1;
    if (q1->n_ext > q2->n_ext)
        return 1;
    if (q1->n_ext < q2->n_ext)
        return -1;
    return 0;
}

static inline void sort_quals(qual_t **quals)
{
    for (int i = 1; i < 4; i++) {
        for (int k = i; k > 0 && cmp_quals(quals[k], quals[k - 1]) > 0; k--) {
            qual_t *tmp = quals[k];
            quals[k] = quals[k - 1];
            quals[k - 1] = tmp;
        }
    }
}

static void compute_qual_freqs(int mer_len, int *reads_lens, char *reads_nts, char *reads_quals, 
                               int nreads, htable_t full_info) 
{
    //int max_read_len = 0;
    int mer_plus = mer_len + 1;
    for (int i = 0; i < nreads; i++) {
        int start = 0;
        int nts_len = reads_lens[i];
        //if (nts_len > max_read_len)
        //    max_read_len = nts_len;
        if (mer_len >= nts_len) 
            continue;

        int subseq_len = nts_len;
        int subseq_pos = 0;
        while (get_valid_base_str(GET_READ(reads_nts, i) + subseq_pos, mer_plus, &start, &subseq_len)) {
            char *nts = strndup(GET_READ(reads_nts, i) + subseq_pos + start, subseq_len);
            char *quals = strndup(GET_READ(reads_quals, i) + subseq_pos + start, subseq_len);

            for (int j = 0; j < subseq_len - mer_len; j++) {
                char *mer = strndup(nts + j, mer_len);
                char extension = (nts + j)[mer_len];

                int offs = j + mer_plus - 1;
                int q = (quals[offs] - _cfg.qual_offset) / 10;
                if (q > 3)
                    q = 3;

                qual_freqs_t *qual_freqs = htable_get_qual_freqs(full_info, mer);
                if (!qual_freqs) {
                    qual_freqs = calloc_chk(1, sizeof(qual_freqs_t));
                    qual_freqs->mer = mer;
                    CHECK_ERR(htable_put_qual_freqs(full_info, mer, qual_freqs, NO_CHECK_DUPS));
                } else {
                    free(mer);
                }
                int bi = get_base_index(extension);
                if (bi < 0) {
                    tprintf("[%d:%d] %s\n", i, nts_len, GET_READ(reads_nts, i));
                    DIE("Invalid base at %d: %d\n", offs, extension);
                }
                qual_freqs->freqs[bi][q]++;
            }
            free(nts);
            free(quals);

            subseq_pos += (start + subseq_len);
            if (subseq_pos >= nts_len) 
                break;
            subseq_len = nts_len - subseq_pos;
        }
    }
}

static void analyse_qual_freqs(htable_t full_info, htable_t mers) 
{
    //Full analysis of quality/frequency profile
    qual_freqs_t *qual_freqs;
    htable_iter_t iter = htable_get_iter(full_info);
    while ((qual_freqs = htable_get_next_qual_freqs(full_info, iter)) != NULL) {
        qual_t *quals[4];
        int n_total = 0;
        for (int i = 0; i < 4; i++) {
            quals[i] = malloc_chk(sizeof(qual_t));
            quals[i]->base = BASES[i];
            // ignore q<10 bases
            quals[i]->n_ext = qual_freqs->freqs[i][1] + qual_freqs->freqs[i][2] + 
                qual_freqs->freqs[i][3];
            quals[i]->n_hi_q = qual_freqs->freqs[i][2] + qual_freqs->freqs[i][3];
            n_total += quals[i]->n_ext;
            quals[i]->rating = categorize_extension(qual_freqs->freqs[i]);
        }

        sort_quals(quals);

        //Rules for choosing next base
        //ratings:
        //0 = No votes
        //1 = One vote
        //2 = nVotes < minViable
        //3 = minDepth > nVotes >= minViable, nHiQ < minViable
        //4 = minDepth > nVotes >= minViable ; nHiQ >= minViable
        //5 = nVotes >= minDepth ; nHiQ < minViable
        //6 = nVotes >= minDepth ; minViable <= nHiQ < minDepth 
        //7 = nHiQ >= minDepth 

        int top_rating = quals[0]->rating;
        int runner_up = quals[1]->rating;
        int top_rated_base = quals[0]->base;

        mer_base_t *mer_base = calloc_chk(1, sizeof(mer_base_t));
        mer_base->mer = strdup(qual_freqs->mer);
        if (top_rating < 3) {         // must have at least minViable bases
            //mer_base->base = 'X';
        } else if (top_rating == 3) {    // must be uncontested   
            if (runner_up == 0) 
                mer_base->base = top_rated_base;
            else 
                mer_base->base = 0;
        } else if (top_rating < 6) {
            if (runner_up < 3)
                mer_base->base = top_rated_base;
            else 
                mer_base->base = 0;
        } else if (top_rating == 6) {  // viable and fair hiQ support
            if (runner_up < 4) 
                mer_base->base = top_rated_base;
            else 
                mer_base->base = 0;
        } else {                     // strongest rating trumps
            if (runner_up < 7) {       
                mer_base->base = top_rated_base;
            } else {
                int k = 0;
                for (int b = 0; b < 4; b++) {
                    if (quals[b]->rating == 7) {
                        mer_base->fork_bases[k] = quals[b]->base;
                        mer_base->fork_n_ext[k] = quals[b]->n_ext;
                        k++;
                    } else {
                        break;
                    }
                }
            }
        }
        CHECK_ERR(htable_put_mer_base(mers, mer_base->mer, mer_base, CHECK_DUPS));
        for (int i = 0; i < 4; i++)
            free(quals[i]);
    }
    free(iter);
}

static char *get_fork_str(char *s, mer_base_t *b) 
{
    strcpy(s, "F");
    for (int i = 0; i < 4; i++) {
        if (!b->fork_bases[i])
            break;
        char buf[100];
        sprintf(buf, "%c%d", b->fork_bases[i], b->fork_n_ext[i]);
        strcat(s, buf);
    }
    return s;
}

#define EXTEND_WALK(walk, base, len)                                \
    do {                                                            \
        if (len >= MAX_WALK_LEN) DIE("walk too long, %ld\n", len); \
        strncat(walk, &(base), 1);                                  \
    } while (0)

static int walk_mers(char *primer1, char *walk_result, htable_t mers, htable_t full_info,
                     char *walk, char *true_primer2) 
{
    char *step_buf = malloc_chk(MAX_WALK_LEN);
    char *step = step_buf;
    strcpy(step, primer1);
    long step_len = strlen(step);
    int success = 0;
    int n_forks = 0;
    int n_steps = 0;
    htable_t loop_check = create_htable(3500, "loop_check");
    int val = 1;

    while (1) {
        if (htable_get_int32(loop_check, step)) {
            strcpy(walk_result, "R");
            break;
        } else {
            CHECK_ERR(htable_put_int32(loop_check, strdup(step), &val, NO_CHECK_DUPS));
        }
        mer_base_t *next = htable_get_mer_base(mers, step);
        if (next) {
            n_steps++;
            if (_cfg.verbose > 2) {
                tprintf("%d : %s->", n_steps, step);
                if (next->base) {
                    tprintf("%c\t", next->base);
                } else if (next->fork_bases[0]) {
                    char buf[100];
                    tprintf("%s\t", get_fork_str(buf, next));
                } else {
                    tprintf("X\t");
                }
                qual_freqs_t *qf = htable_get_qual_freqs(full_info, step);
                for (int q = 0; q < 4; q++) 
                    tprintf("[%d %d %d %d]", qf->freqs[q][0], qf->freqs[q][1], qf->freqs[q][2], 
                            qf->freqs[q][3]);
                tprintf("\n");
            }
            if (next->base) {
                step++;
                step_len++;
                EXTEND_WALK(step, next->base, step_len);
                EXTEND_WALK(walk, next->base, strlen(walk));
                if (endswith(walk, true_primer2) == 0) {
                    success = 1;
                    break;
                }
            } else if (_cfg.poly_mode && !n_forks && next->fork_bases[0]) {
                //  Biallelic positions only (for now) maximum vote path is taken
                char c = 0;
                if (next->fork_bases[1] && !next->fork_bases[2]) {
                    c = next->fork_n_ext[0] > next->fork_n_ext[1] ? 
                        next->fork_bases[0] : next->fork_bases[1];
					tprintf_vvv("Polymorphic conditions met .. "
								"attempting max frequency resolution.\n");
                } else {
                    char buf[100];
                    strcpy(walk_result, get_fork_str(buf, next));
                    break;
                }
                step++;
                step_len++;
                EXTEND_WALK(walk, c, strlen(walk));
                EXTEND_WALK(step, c, step_len);
                n_forks++;
                if (endswith(walk, true_primer2) == 0) {
                    success = 1;
                    walk_result[0] = 0;
                    break;
                }
            } else {
                if (next->fork_bases[0]) {
                    char buf[100];
                    strcpy(walk_result, get_fork_str(buf, next));
                } else { 
                    strcpy(walk_result, "X");
                }
                break;
            }
        } else {
            strcpy(walk_result, "X");
            break;
        }
    }
    free(step_buf);
    destroy_htable(loop_check, FREE_KEYS);
    return success;
}

static void bridge_iter(char *c1_seq, char *c2_seq, int *seq_ref_lens, char *seq_ref_nts, 
                        char *seq_ref_quals, int nreads, int dirn, char **dwalk, char *dfail)
{
	tprintf_vvv("Attempting %s bridge: (%d read(s) available)\n", DIRNS[dirn], nreads);

    int *reads_lens = malloc_chk(sizeof(int) * nreads);
    char *reads_nts = malloc_chk((_cfg.max_read_len + 1) * nreads);
    char *reads_quals = malloc_chk((_cfg.max_read_len + 1) * nreads);
    int max_read_length = 0;

    for (int i = 0; i < nreads; i++) {
        reads_lens[i] = seq_ref_lens[i];
        strcpy(GET_READ(reads_nts, i), GET_READ(seq_ref_nts, i));
        if (reads_lens[i] > max_read_length)
            max_read_length = reads_lens[i];
        strcpy(GET_READ(reads_quals, i), GET_READ(seq_ref_quals, i));
        if (dirn == DIRN_LEFT) {
            switch_code(reverse(GET_READ(reads_nts, i)));
            reverse(GET_READ(reads_quals, i));
        }
    }

    const int min_mer_len = 13;
    int mer_len = _cfg.mer_size;
    int downshift = 0;
    int upshift = 0;

    char *max_walk = NULL;
    int max_walk_len = 0;

    char *c1s, *c2s;
    if (dirn == DIRN_LEFT) {
        c1s = c2_seq;
        c2s = c1_seq;
        switch_code(reverse(c1s));
        switch_code(reverse(c2s));
    } else {
        c1s = c1_seq;
        c2s = c2_seq;
    }

    int sl1 = strlen(c1s);

    char *true_primer2 = strndup(c2s, _cfg.mer_size);
    char *primer1 = NULL;
    char *primer2 = NULL;
    char walk_result[100] = "";
    while (1) {
        walk_result[0] = 0;
        // Allow k-mers less than mer_size to be used 
        if (mer_len < _cfg.mer_size) {
            int the_rest = _cfg.mer_size - mer_len;
            primer1 = strndup(c1s + sl1 - _cfg.mer_size, mer_len);
            primer2 = strndup(c2s + the_rest, mer_len);
        } else {
            primer1 = strdup(c1s + sl1 - mer_len);
            primer2 = strndup(c2s, mer_len);
        }
        if (strlen(primer1) < mer_len) {
            free(primer1);
            primer1 = NULL;
        }
        if (strlen(primer2) < mer_len) {
            free(primer2);
            primer2 = NULL;
        }
        if (!primer1 || !primer2) {
			tprintf_vvv("Unable to find k=%d seeds\n", mer_len);
            if (!max_walk) {
                max_walk = malloc_chk(1);
                max_walk[0] = 0;
                strcpy(walk_result, "X");
                strcpy(dfail, walk_result);
            }
            break;
        }

        htable_t full_info = create_htable(_max_htable_entries, "full_info");
        compute_qual_freqs(mer_len, reads_lens, reads_nts, reads_quals, nreads, full_info);
		int nentries = htable_num_entries(full_info);
		if (nentries > htable_capacity(full_info) * 0.65) {
			_max_htable_entries = 2 * nentries;
			tprintf("Expanding max htable entries to %d\n", _max_htable_entries);
		}
        htable_t mers = create_htable(_max_htable_entries, "mers");
        analyse_qual_freqs(full_info, mers);

        char *walk = malloc_chk(MAX_WALK_LEN);
        strcpy(walk, primer1);

        if (!max_walk) {
            max_walk = strdup(walk);
            max_walk_len = strlen(max_walk);
            strcpy(dfail, walk_result);
        }

        int success = walk_mers(primer1, walk_result, mers, full_info, walk, true_primer2);

        // clean up
        destroy_htable(full_info, FREE_KEYS|FREE_VALS);
        destroy_htable(mers, FREE_KEYS|FREE_VALS);

        // Trim off extra lead bases if upshifted
        int additional_bases = mer_len - _cfg.mer_size;
        if (additional_bases < 0) 
            additional_bases = 0;
        int walk_len = strlen(walk) - additional_bases;
        if (success || (walk_len > max_walk_len)) {
            if (max_walk)
                free(max_walk);
            max_walk = strdup(walk + additional_bases);
            max_walk_len = walk_len;
            strcpy(dfail, walk_result);
        }

        if (primer1)
            free(primer1);
        if (primer2)
            free(primer2);
        free(walk);
        primer1 = NULL;
        primer2 = NULL;

        if (walk_result[0] == 'F' || walk_result[0] == 'R') {
            mer_len += 2;
            upshift = 1;
            if (downshift == 1 || mer_len >= max_read_length) 
                break;
			tprintf_vvv("Degeneracy encountered; upshifting (k->%d)\n", mer_len);
        } else if (walk_result[0] == 'X') {
            mer_len -= 2;
            downshift = 1;
            if (upshift == 1 || mer_len < min_mer_len) 
                break;
			tprintf_vvv("Termination encountered; downshifting (k->%d)\n", mer_len);
        } else {
            break;
        }
    }

    *dwalk = max_walk;

    tprintf_vv("MAX%s [%s%s]\n", dirn == DIRN_RIGHT ? "RIGHT" : "LEFT", 
               *dwalk ? *dwalk : "", dfail);

    free(true_primer2);
    free(reads_lens);
    free(reads_nts);
    free(reads_quals);
    return;
}

static char *patch(char *right_walk, char *left_walk, int gap_size, double gap_uncertainty) 
{
    int min_acceptable_overlap = 10;
    int right_len = strlen(right_walk);
    int left_len = strlen(left_walk);
    int gu = bankers_round(gap_uncertainty);
    int ideal_len = 2 * _cfg.mer_size + gap_size;
    int chop_left = left_len - ideal_len;
    int chop_right = right_len - ideal_len;

	tprintf_vvv("Attempting patch [%d][%d] (gap: %d +/- %.1f)\n", 
                right_len, left_len, gap_size, gap_uncertainty);

    char *rwalk = strdup(right_walk);
    char *lwalk = strdup(left_walk);
    char *test_right = rwalk;
    char *test_left = lwalk;
    if (chop_left > 0) 
        test_left += chop_left;
    int tl_len = strlen(test_left);
    if (chop_right > 0) 
        test_right[right_len - chop_right] = 0;
    int tr_len = strlen(test_right);

    int max_buf = MAX_BEST_GUESS_SEQ_LEN;
    char *best_guess_seq = malloc_chk(max_buf);
    best_guess_seq[0] = 0;
    double best_guess_delta = gap_uncertainty + 1;
    int n_best_guesses = 0;
    int n_guesses = 0;
    for (int o = -gu; o <= gu; o++) {
        int overlap = tr_len + tl_len - (ideal_len + o);
        if (overlap < min_acceptable_overlap || overlap > tr_len || overlap > tl_len) 
            continue;
        char *p1_suffix = test_right + tr_len - overlap;
        char *p2_suffix = test_left + overlap;
        if (strncmp(p1_suffix, test_left, overlap) != 0) 
            continue;
        int delta = abs(o);
        n_guesses++;
        int newlen = strlen(test_right) + strlen(p2_suffix);
        if (newlen >= max_buf) {
            WARN("buf for best guess seq len %d is too small\n", newlen);
            max_buf = newlen + 2;
            best_guess_seq = realloc_chk(best_guess_seq, max_buf);
        }
        if (delta < best_guess_delta) {
            sprintf(best_guess_seq, "%s%s", test_right, p2_suffix);
            best_guess_delta = delta;
            n_best_guesses = 1;
        } else if (delta == best_guess_delta) {
            if (strlen(test_right) + strlen(p2_suffix) < strlen(best_guess_seq)) {
                sprintf(best_guess_seq, "%s%s", test_right, p2_suffix);
                best_guess_delta = delta;
            }
            n_best_guesses++;
        }
    }
    free(rwalk);
    free(lwalk);

    if (n_guesses) {
        tprintf_vv("%d potential patches identified.  "
				   "%d best guess(es) differ from gap estimate by %.0f\n",
				   n_guesses, n_best_guesses, best_guess_delta);
        return best_guess_seq;
    }
    free(best_guess_seq);
    tprintf_vv("No valid patches found.\n");
    return NULL;
}

static char *mer_walk_dirn(gap_t *gap, char *c1_seq, char *c2_seq, Buffer report_note, char **dwalk, 
                           char *dfail, int *check, int dirn) 
{
    _nwalks++;
    START_TIMER(T_MER_WALKS);
    *dwalk = NULL;
    dfail[0] = 0;
    char *closure = NULL;
    bridge_iter(c1_seq, c2_seq, (int*)gap->reads_lens, (char*)gap->reads_nts, (char*)gap->reads_quals, 
                gap->nreads, dirn, dwalk, dfail);
    if (strcmp(*dwalk, "") == 0)
        DIE("empty dwalk\n");
    if (!*dwalk) 
        return NULL;
    if (dirn == DIRN_LEFT) {
        switch_code(reverse(*dwalk));
        if (dfail[0] == 'F')
            switch_code(dfail);
    }
    if (!dfail[0]) {
        *check = check_closure(*dwalk, gap->primer1, gap->primer2, gap->size, gap->uncertainty);
        add_to_report(report_note, "%sClosure=%s;%sCheck=%d;", DIRNS[dirn], *dwalk, DIRNS[dirn], 
                      *check);
        if (*check == 0)
            closure = strdup(*dwalk);
    }
    stop_timer(T_MER_WALKS);
    return closure;
}

static void lowercase_gap(gap_t *gap, char *closure, int clen) 
{
    // For more clear output - lowercase the sequence that was closed
    //  Lower case as little of the primer sequences as possible:
    int min_gap_mask = 2 * 5;
    int p1len = strlen(gap->primer1);
    int p2len = strlen(gap->primer2);
    int start_gap = p1len;
    int end_gap = clen - p2len;
    if (clen < p1len + p2len + min_gap_mask) {
        start_gap = clen / 2 - min_gap_mask / 2;
        end_gap = clen / 2 + min_gap_mask / 2;
        if (clen % 2)
            end_gap++;
    }
    for (int i = start_gap; i < end_gap; i++) 
        closure[i] = tolower(closure[i]);
}

static void close_gaps(gap_t *gaps, int ngaps)
{
    START_TIMER(T_CLOSE_GAPS);
    serial_printf("Closing gaps... ");
    tprintf_flush("Attempting to close %d gaps...\n", ngaps);

    int success = 0;
	int n_span_closures = 0;
	int n_left_closures = 0;
	int n_right_closures = 0;
	int n_patch_closures = 0;
    char buf[MAX_BUF];
    Buffer report_note = initBuffer(MAX_BUF*2), report_line = initBuffer(MAX_BUF);
#ifdef CONFIG_SHOW_PROGRESS
    int tick_size = (ngaps+10-1) / 10;
#endif
    for (int i = 0; i < ngaps; i++) {
        // for debugging
#ifdef SOLO_GAP
        if (i != SOLO_GAP) 
            continue;
#endif
#ifdef CONFIG_SHOW_PROGRESS
        if (i && !(i % tick_size)) {
            serial_printf("%d ", i / tick_size);
            tprintf_flush("%d ", i / tick_size);
        }
#endif
        gap_t *gap = (gap_t*)&gaps[i];
        gap->closure = NULL;
        if (!gap->nreads) 
            continue;
        resetBuffer(report_note);
        shared contig_t *contig1;
        get_contig_chk(contig1, gap->contig1_id);
        shared contig_t *contig2;
        get_contig_chk(contig2, gap->contig2_id);
        shared scaffold_t *scaffold = get_scaffold(gap->scaff_id);
        if (!contig1 || !contig2 || !scaffold) {
            DIE("Missing for gap %d: !contig1 %d !contig2 %d !scaffold %d, skipping...\n", 
                 gap->id, !contig1, !contig2, !scaffold);
            //continue;
        }
        if (!contig1->seq || !contig2->seq) {
            DIE("Missing contigs for gap %d: !contig1 %d !contig2 %d, skipping...\n", 
                 gap->id, !contig1->seq, !contig2->seq);
            continue;
        }

        // compute the primers here
        gap->primer1 = get_primer(contig1, gap->contig_ext1, 5);
        gap->primer2 = get_primer(contig2, gap->contig_ext2, 3);
        resetBuffer(report_line);
        get_report_line(gap, report_line);
        if (_cfg.exclude_repeats) {
            if (scaffold->depth > _cfg.exclude_repeats) {
                tprintf_vv("\n*************\nRepeat scaffold gap excluded. (depth = %.6f) %s\n",
						   scaffold->depth, get_gap_info(gap, buf, MAX_BUF - 1));
                tprintf_vv("%s\tFAILED\tscaffDepth=%.6f", report_line->buf, scaffold->depth);
                continue;
            }
        }

        tprintf_vv("\n*************\nAttempt to close %d: %s\n",  
				   i, get_gap_info(gap, buf, MAX_BUF - 1));

        gap->uncertainty = 3 * gap->uncertainty + 5.5;

        int span_check;
        char *closure = splinting_reads(gap, &span_check, report_note);
        // If splints fail try a mer-walk
        int right_check = 0;
        char *right_walk = NULL;
        char right_fail[100] = "";
        int left_check = 0;
        char *left_walk = NULL;
        char left_fail[100] = "";

        if (closure) {
			n_span_closures++;
        } else {
            char *c1_seq = get_seq(contig1, gap->primer1, DIRN_LEFT);
            char *c2_seq = get_seq(contig2, gap->primer2, DIRN_RIGHT);
            if (c1_seq && c2_seq) {
                closure = mer_walk_dirn(gap, c1_seq, c2_seq, report_note, &right_walk, right_fail, 
                                        &right_check, DIRN_RIGHT);
                if (closure) {
					n_right_closures++;
				} else {
                    closure = mer_walk_dirn(gap, c1_seq, c2_seq, report_note, &left_walk, left_fail,
                                            &left_check, DIRN_LEFT);
					if (closure)
						n_left_closures++;
				}
            }
            if (c1_seq)
                free(c1_seq);
            if (c2_seq)
                free(c2_seq);
		}
        snprintf(buf, MAX_BUF - 1, "[%d\t%s\t%d\t%s]", gap->nreads, 
                 gap->primer1, gap->size, gap->primer2);

        START_TIMER(T_PATCHING);
        // If walks failed to close, try to patch between left and right walk
        if (!closure && right_fail[0] && left_fail[0]) {
            closure = patch(right_walk, left_walk, gap->size, gap->uncertainty);
            if (closure) {
                int patch_check = check_closure(closure, gap->primer1, gap->primer2, gap->size, 
                                                gap->uncertainty);
                add_to_report(report_note, "patchClosure=%s;patchCheck=%d;", closure, patch_check);
                if (patch_check != 0) {
                    free(closure);
                    closure = NULL;
                } else {
					n_patch_closures++;
				}
            }
        }
        stop_timer(T_PATCHING);

        if (closure) {
            int closed_gap_size = strlen(closure) - 2 * _cfg.mer_size;
            success = 1;
            _nsuccess++;
            tprintf_vv("%s successfully closed gap %d: %d %s\n", 
					   buf, i, closed_gap_size, closure);
#ifdef CONFIG_SANITY_CHECK            
            if (strncmp(closure, gap->primer1, strlen(gap->primer1)) != 0)
                DIE("closure for gap %d doesn't start with primer1\n", i);
            if (endswith(closure, gap->primer2) != 0)
                DIE("closure for gap %d doesn't end with primer2\n", i);
#endif
            lowercase_gap(gap, closure, strlen(closure));
            int len_closure = strlen(closure);

            if (len_closure < 4)
                WARN("short closure len %d: %s\n", len_closure, closure);

            gap->closure = upc_alloc(len_closure + 1);
            upc_memput(gap->closure, closure, len_closure + 1);
            gap->len_closure = len_closure;
        } else {
            success = 0;
            _nfailure++;
            tprintf_vv("%s failed to close gap %d (%d;%d:%s;%d:%s)\n",
					   buf, i, span_check, right_check, right_fail, left_check, left_fail);
        }

        if (right_walk)
            free(right_walk);
        if (left_walk)
            free(left_walk);
        if (closure)
            free(closure);

        tprintf_vv("%s\n", report_line->buf);
        tprintf_vv("%s\t%s\n", success ? "SUCCESS" : "FAILED", report_note->buf);
        tprintf_v("%d: %s\t%s\n", i, success ? "SUCCESS" : "FAILED", report_note->buf);
    }

    stop_timer(T_CLOSE_GAPS);
    freeBuffer(report_line);
    freeBuffer(report_note);
    tprintf_flush("...  successfully closed %d gaps (%d failed to close)\n", 
                  _nsuccess, _nfailure);
	tprintf_flush("Closures: %d spans, %d right, %d left, %d patched\n", 
				  n_span_closures, n_right_closures, n_left_closures, n_patch_closures);
    tprintf_flush("Number of walks: %d\n", _nwalks);
    //double t_all = get_timer_all_max(T_CLOSE_GAPS);
    //serial_printf("done in %.3f s\n",  t_all);
    tprintf_flush("done in %.3f s\n", get_elapsed_time(T_CLOSE_GAPS));
    UPC_TIMED_BARRIER;
    serial_printf("done in %.3f s\n", get_elapsed_time(T_CLOSE_GAPS));
}

static void print_closures(gap_t *gaps, int ngaps) 
{
    char closures_fname[MAX_FILE_PATH];
    serial_printf("printing closures...");
    sprintf(closures_fname, "%s/closures.%d.out", _cfg.base_dir, MYTHREAD);
    get_rank_path(closures_fname, MYTHREAD);
    FILE *closures_file = fopen_chk(closures_fname, "w");
    for (int i = 0; i < ngaps; i++) {
        gap_t *gap = (gap_t*)&gaps[i];
        if (gap->closure) {
            fprintf(closures_file, "Scaffold%d\tContig%d.%d\t%s\tContig%d.%d\t\%s\t%s\n",
                    gap->scaff_id, gap->contig1_id, gap->contig_ext1, gap->primer1, 
                    gap->contig2_id, gap->contig_ext2, gap->primer2, (char*)gap->closure);
        }
    }            
    fclose(closures_file);
    serial_printf("done.\n");
}


static void print_fasta(FILE *f, int scaff_id, char *seq, int seq_len)
{
    fprintf(f, ">Scaffold%d\n", scaff_id);
    int nlines = INT_CEIL(seq_len, FASTA_LINE_WIDTH);
    char buf[FASTA_LINE_WIDTH + 1];
    for (int i = 0; i < nlines; i++) {
        strncpy(buf, seq + (i * FASTA_LINE_WIDTH), FASTA_LINE_WIDTH);
        buf[FASTA_LINE_WIDTH] = 0;
#ifdef CHECK_CHARS_IN_OUTPUT
        CHECK_SEQ_POS(buf);
#endif
        fprintf(f, "%s\n", buf);
    }
}

static char *to_nearest_units(char *s, long val)
{
    const double one_gb = 1000000000;
    const double one_mb = 1000000;
    const double one_kb = 1000;
    if (val >= one_gb)
        sprintf(s, "%.1f GB", (double)val / one_gb);
    else if (val >= one_mb)
        sprintf(s, "%.1f MB", (double)val / one_mb);
    else if (val >= one_kb)
        sprintf(s, "%.1f KB", (double)val / one_kb);
    else 
        sprintf(s, "%ld", val);
    return s;
}

int cmp_int(const void *x, const void *y)
{
    int xi = *(int*)x;
    int yi = *(int*)y;
    if (xi < yi) return 1;
    if (xi > yi) return -1;
    return 0;
}

static void calc_n50(FILE *f, shared int *lens, int n, int max_i, const char *desc)
{
    if (!n)
        SDIE("Cannot calculate assembly stats for 0 elements\n");
    UPC_TIMED_BARRIER;
    // just sort in thread 0
    if (MYTHREAD == 0) {
        long tot_size = 0;
        int *local_lens = malloc_chk(sizeof(int) * n);
        int j = 0;
        for (int i = 0; i < max_i; i++) {
            if (!lens[i]) 
                continue;
            local_lens[j] = lens[i];
            tot_size += local_lens[j];
            j++;
        }
        // now use system qsort
        qsort(local_lens, n, sizeof(int), cmp_int);
        //for (int i = 0; i < n; i++)
        //    dbg("%d\t%d\n", i, local_lens[i]);
        // now run the algorithm
        long running_tot = 0;
        int n_index = 1;
        int N50 = 0;
        int L50 = 0;
        char buf[100];
        for (int i = 0; i < n; i++) {
            running_tot += local_lens[i];
            while (running_tot > tot_size * n_index / 100) {
                if (n_index == 50) {
                    L50 = local_lens[i];
                    N50 = i + 1;
                }
                n_index++;
            }
        }
        free(local_lens);
        if (!MYTHREAD)
            fprintf(f, "\tMain genome %s N/L50: %d/%s\n", desc, N50, to_nearest_units(buf, L50));
    }
}

static void calc_assembly_stats(assembly_stats_t *stats)
{
    START_TIMER(T_CALC_ASSEMBLY_STATS);
    FILE *f = fopen_chk("final_assembly.stats", "w");
    char buf[100];
    if (!MYTHREAD)
        fprintf(f, "Assembly stats (excluding scaffolds < %d):\n", _cfg.scaff_len_cutoff);
    stats->scaff_tot = reduce_int(stats->scaff_tot, UPC_ADD);
    if (!MYTHREAD)
        fprintf(f, "\tMain genome scaffold total: %s\n", to_nearest_units(buf, stats->scaff_tot));
    stats->contig_tot = reduce_int(stats->contig_tot, UPC_ADD);
    if (!MYTHREAD)
        fprintf(f, "\tMain genome contig total:   %s\n", to_nearest_units(buf, stats->contig_tot));
    stats->scaff_len_tot = reduce_long(stats->scaff_len_tot, UPC_ADD);
    if (!MYTHREAD)
        fprintf(f, "\tMain genome scaffold sequence total: %s\n", 
                  to_nearest_units(buf, stats->scaff_len_tot));
    stats->contig_len_tot = reduce_long(stats->contig_len_tot, UPC_ADD);
    if (!MYTHREAD)
        fprintf(f, "\tMain genome contig sequence total:   %s (->  %.1f%% gap)\n", 
                  to_nearest_units(buf, stats->contig_len_tot),
                  stats->scaff_len_tot == 0 ? 0.0 : 
                  100.0 * (1.0 - (double)stats->contig_len_tot / stats->scaff_len_tot));
    calc_n50(f, stats->scaff_lens, stats->scaff_tot, get_max_scaffolds(), "scaffold");
    //calc_n50(f, stats->contig_lens, stats->contig_tot, get_max_contigs(), "contig");

    fclose(f);
    stop_timer(T_CALC_ASSEMBLY_STATS);
    serial_printf("done in %.4f s\n", get_elapsed_time(T_CALC_ASSEMBLY_STATS));
}

static void check_scaff_seq_len(char **scaff_seq, int len, int *max_scaff_seq_len)
{
    if (len < 0)
        DIE("Scaffold sequence < 0, i.e. %d\n", len);
    if (len > *max_scaff_seq_len) {
        *max_scaff_seq_len = 2 * len;
        *scaff_seq = realloc_chk(*scaff_seq, *max_scaff_seq_len);
    }
}

static void scaff_to_fasta(assembly_stats_t *stats, shared gap_set_t *gaps, int tot_gaps)
{
    tprintf_flush("Converting scaffolds and gaps to fasta...");
    serial_printf("Converting scaffolds and gaps to fasta...");
    START_TIMER(T_SCAFF_TO_FASTA);
    char fname[MAX_FILE_PATH];
    sprintf(fname, "%s/assembly.%d.fa", _cfg.base_dir, MYTHREAD);
    get_rank_path(fname, MYTHREAD);
    FILE *f = fopen_chk(fname, "w");
    FILE *f_short = NULL;
    if (_cfg.scaff_len_cutoff) {
        sprintf(fname, "%s/assembly-%d.%d.fa", _cfg.base_dir, _cfg.scaff_len_cutoff, MYTHREAD);
        get_rank_path(fname, MYTHREAD);
        f_short = fopen_chk(fname, "w");
    }

    shared scaffold_t *scaffolds = get_scaffolds();
    int max_scaffolds = get_max_scaffolds();

    // statistics
    stats->scaff_lens = upc_all_alloc(max_scaffolds + 1, sizeof(int));

    char *contig_seq = malloc_chk(MAX_SCAFF_SEQ_LEN);
    char *closure = malloc_chk(MAX_SCAFF_SEQ_LEN);
    strcpy(closure, "");
    int max_scaff_seq_len = MAX_SCAFF_SEQ_LEN;
    char *scaff_seq = malloc_chk(max_scaff_seq_len);

    int slen, clen, offset;
    shared scaffold_t *scaffold;
    shared contig_t *c1, *c2;
    int s_start, gap_size, gap_excess, cutback1, cutback2, ncount, breaks, scaff_seq_len;
    upc_forall (int i = 0; i < max_scaffolds; i++; i) {
        scaff_seq[0] = 0;
        slen = 0;
        clen = 0;
        offset = 0;
        scaffold = &scaffolds[i];
        if (scaffold->id == -1)
            continue;
        for (int j = 0; j < scaffold->ncontigs; j++) {
            get_contig_chk(c1, scaffold->contig_ids[j]);
#ifdef DBG_SCAFF_FASTA            
            dbg("Scaffold%d\tContig%d %c\n", scaffold->id, c1->id, c1->strand);
#endif
            UPC_MEMGET_STR(contig_seq, c1->seq, c1->seq_len + 1);
            if (c1->strand == '-') 
                switch_code(reverse(contig_seq));
            
            clen = strlen(contig_seq) - offset;
            check_scaff_seq_len(&scaff_seq, clen + slen, &max_scaff_seq_len);
            strcat(scaff_seq, contig_seq + offset);
            slen += clen;
            offset = 0;
            
            shared [] gap_t *gap = NULL;
            if (c1->strand == '-' && c1->gapi_3 != -1) 
                gap = get_gap(c1->gapi_3, __LINE__, __FILE__, ABORT_ON_FAIL);
            else if (c1->strand == '+' && c1->gapi_5 != -1) 
                gap = get_gap(c1->gapi_5, __LINE__, __FILE__, ABORT_ON_FAIL);

            int gap_closure = 0;
            if (gap != NULL) {
                if (gap->closure != NULL) {
                    UPC_MEMGET_STR(closure, gap->closure, gap->len_closure + 1);
#ifdef DBG_SCAFF_FASTA            
                    *c2 = NULL;
                    if (gap->contig1_id == c1->id) 
                        get_contig_chk(c2, gap->contig2_id);
                    else
                        get_contig_chk(c2, gap->contig1_id);
                    dbg("    closure Contig%d:Contig%d\n", c1->id, c2->id); 
                    dbg("    closure %s, %d\n", closure, strlen(closure));
#endif
                    clen = strlen(closure);
                    check_scaff_seq_len(&scaff_seq, slen - _cfg.mer_size + clen, &max_scaff_seq_len);
                    strcpy(scaff_seq + slen - _cfg.mer_size, closure);
                    slen += (clen - _cfg.mer_size);
                    offset = _cfg.mer_size;
                    gap_closure = 1;
                }
            }
            if (!gap_closure) {
                // fill in with NNs 
                if (j + 1 < scaffold->ncontigs) {
                    get_contig_chk(c2, scaffold->contig_ids[j + 1]);
                    s_start = c2->s_start;
                    gap_size = s_start - c1->s_end - 1;
#ifdef DBG_SCAFF_FASTA            
                    dbg("    seq no closure, gap size %d, start %d, coord %d\n", 
                        gap_size, s_start, c1->s_end);
#endif
                    // Trim contigs back around residual negative gaps to eliminate potential 
                    // redundancy
                    if (gap_size < _cfg.min_gap_size) {
                        gap_excess = _cfg.min_gap_size - gap_size;
                        cutback1 = !(gap_excess % 2) ? gap_excess / 2 : (gap_excess + 1) / 2;
                        cutback2 = !(gap_excess % 2) ? gap_excess / 2 : (gap_excess - 1) / 2;
#ifdef DBG_SCAFF_FASTA            
                        dbg("    Nibbling back around gap: Contig%d:Contig%d\n", c1->id, c2->id);
#endif
                        slen -= cutback1;
                        check_scaff_seq_len(&scaff_seq, slen, &max_scaff_seq_len);
                        scaff_seq[slen] = 0;
                        if (cutback2 + _cfg.mer_size >= c2->seq_len && slen >= cutback2) {
                            // Don't trim the upcoming contig if the trimming would leave less than
                            // a primer-size of sequence.  We will need at least that much sequence
                            // to prime the next gap.  Instead, cut back more from the already
                            // processed scaffold
                            slen -= cutback2;
                            check_scaff_seq_len(&scaff_seq, slen, &max_scaff_seq_len);
                            scaff_seq[slen] = 0;
                        } else {
                            offset = cutback2;
                        }
                        gap_size = _cfg.min_gap_size;
                    }
                    check_scaff_seq_len(&scaff_seq, slen + gap_size, &max_scaff_seq_len);
                    memset(scaff_seq + slen, 'N', gap_size);
                    slen += gap_size;
                    scaff_seq[slen] = 0;
                }
            }
        }
        scaff_seq_len = strlen(scaff_seq);
        if (scaff_seq_len >= _cfg.scaff_len_cutoff) {
            print_fasta(f, scaffold->id, scaff_seq, scaff_seq_len);
            // gather assembly stats
            stats->scaff_lens[scaffold->id] = scaff_seq_len;
            stats->scaff_tot++;
            // compute number of contigs by splitting over at least 25 consecutive N's
            ncount = 0;
            breaks = 0;
            for (int j = 0; j < scaff_seq_len; j++) {
                if (scaff_seq[j] == 'N') {
                    ncount++;
                } else {
                    if (ncount >= 25) {
                        stats->contig_tot++;
                        breaks += ncount;
                    }
                    ncount = 0;
                }
            }
            // the scaffold itself is a contig
            stats->contig_tot++;
            stats->contig_len_tot += (scaff_seq_len - breaks);
            stats->scaff_len_tot += scaff_seq_len;
        } else {
            print_fasta(f_short, scaffold->id, scaff_seq, scaff_seq_len);
            stats->scaff_lens[scaffold->id] = 0;
        }
    }
    fclose(f);
    if (_cfg.scaff_len_cutoff) 
        fclose(f_short);
    free(scaff_seq);
    free(contig_seq);
    free(closure);
    UPC_TIMED_BARRIER;
    stop_timer(T_SCAFF_TO_FASTA);
    tprintf_flush("done in %.4f s\n", get_elapsed_time(T_SCAFF_TO_FASTA));
    serial_printf("done in %.4f s\n", get_elapsed_time(T_SCAFF_TO_FASTA));
}

int main(int argc, char** argv)
{
    if (argc == 1)
        return 0;    

    upc_tick_t start_time = upc_ticks_now();

    init_timers();
    START_TIMER(T_OVERALL);
    get_config(argc, argv, &_cfg);
#ifdef CONFIG_SANITY_CHECK
    tprintf_flush("  sanity checking enabled\n");
#endif
#ifdef CONFIG_CHECK_SEQS
    tprintf_flush("  checking final sequences for illegal characters\n");
#endif

    if (strncmp(_cfg.base_dir, "/dev/shm", 8) == 0) {
        // FIXME: crude hack for removing files from shm generated by previous stages
        // This should be somewhere else outside this code
        char prev_name[MAX_FILE_PATH];
        sprintf(prev_name, "/dev/shm/output_%d_bubble_contigs.fasta", MYTHREAD);
        unlink(prev_name);
        tprintf("unlinking %s\n", prev_name);
        sprintf(prev_name, "/dev/shm/output_%d_contigs.fasta", MYTHREAD);
        unlink(prev_name);
        tprintf("unlinking %s\n", prev_name);
        for (int i = 0; i < _cfg.num_lib_names; i++) {
            sprintf(prev_name, "/dev/shm/%s-spans_%d", _cfg.library_names[i], MYTHREAD);
            unlink(prev_name);
            tprintf("unlinking %s\n", prev_name);
            sprintf(prev_name, "/dev/shm/%s-splints_%d", _cfg.library_names[i], MYTHREAD);
            unlink(prev_name);
            tprintf("unlinking %s\n", prev_name);
        }
    }
    shared gap_set_t *gaps = upc_all_alloc(THREADS, sizeof(gap_set_t));
    init_globals(&_cfg, gaps);
    process_scaffold_report_file();
    process_contigs_file();

    load_balance_gaps();

    for (int i = 0; i < _cfg.num_lib_names; i++) {
        process_meraligner_file(_cfg.library_names[i], _cfg.reverse_complement[i], 
                                _cfg.insert_size[i], _cfg.insert_sigma[i], 
                                _cfg.five_prime_wiggle_room[i], _cfg.three_prime_wiggle_room[i]);
        process_fastq_file(_cfg.library_names[i]);
    }

    //print_timers(0);
    
    if (_cfg.print_gaps)
        print_gaps();
    close_gaps((gap_t*)gaps[MYTHREAD].gaps, gaps[MYTHREAD].ngaps);
    if (_cfg.print_closures)
        print_closures((gap_t*)gaps[MYTHREAD].gaps, gaps[MYTHREAD].ngaps);

    int tot_ngaps_closed = reduce_int(_nsuccess, UPC_ADD);
    int tot_ngaps_unclosed = reduce_int(_nfailure, UPC_ADD);
    int tot_gaps = tot_ngaps_closed + tot_ngaps_unclosed;
    double perc_failed = tot_gaps == 0 ? 0.0 : (double)tot_ngaps_unclosed / tot_gaps * 100.0;
    tprintf_flush("Total gaps closed %d (%d failed to close, %.2f %%)\n", 
                  tot_ngaps_closed, tot_ngaps_unclosed, perc_failed);
    serial_printf("Total gaps closed %d (%d failed to close, %.2f %%)\n", 
                  tot_ngaps_closed, tot_ngaps_unclosed, perc_failed);

    //print_timers(0);
    
    assembly_stats_t stats = {0};
    scaff_to_fasta(&stats, gaps, tot_gaps);
    
    stop_timer(T_OVERALL);
    double t_all = get_timer_all_max(T_OVERALL);
    tprintf_flush("Peak memory usage: %d mb\n", get_max_mem_usage_mb());

    serial_printf("Elapsed time for all threads: %.4f s\n", t_all);
    
    print_timers(0);

    UPC_TIMED_BARRIER;
    serial_printf("Overall time for %s is %.2f s\n", basename(argv[0]), 
                  ((double)upc_ticks_to_ns(upc_ticks_now() - start_time) / 1000000000.0));
    return 0;
}

