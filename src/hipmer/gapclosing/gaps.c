/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <upc.h>
#include <upc_nb.h>

#include "utils.h"
#include "gaps.h"
#include "timers.h"
#include "../fqreader/fq_reader.h"
#include "tracing.h"
#include "../common/common.h"
#include "../common/Buffer.h"
#include "../common/upc_compatibility.h"
#include "../common/kseq.h"
#include "dhtable.h"

//#define DBG_GET_ALIGNMENT
//#define DBG_PAIR_ALNS
//#define DBG_PUT_READ_ORIENT
//#define DBG_READ_DATA
//#define DBG_CONTIGS
//#define DBG_SCAFFOLDS
//#define DBG_SEQUENCES
//#define DBG_DREAD_AFTER
//#define DBG_FASTQ

#define SHOW_MEMORY
//#define CHECK_READ_DATA_HASH

// primes for the hash table
// 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
// 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741

#define MAX_READ_GAP_IDS 100
#define MAX_READS_PER_GAP 5000
#define MAX_SEQ_LEN 5000000
#define MAX_BM_LINE_LEN 256
#define MAX_PAIR_ALNS 6
#define STEAL_BLOCK 100

// to create the hash table for reads, multiplied by the number of reads - gives approx load factor
#define READ_HASH_TABLE_LOAD_FACTOR 0.2

#define REJECT_FORMAT 0
#define REJECT_MINLEN 1
#define REJECT_UNINFORMATIVE 2
#define REJECT_5_TRUNCATED 3
#define REJECT_3_TRUNCATED 4
#define REJECT_SINGLETON 5
#define REJECT_NO_ASSOCIATED_GAP 6
#define REJECT_TRUNC5_FAIL 7
#define REJECT_TRUNC3_FAIL 8
#define REJECT_NO_SCAFFOLD_INFO 9
#define REJECT_POTENTIAL_INNIE 10
#define REJECT_POTENTIAL_SHORTY 11
#define REJECT_HEADER 12
#define NUM_REJECT_REASONS (REJECT_HEADER + 1)
#define GOOD_ALIGNMENT 100

#define GET_SCAFF_ID(scaff) (atoi((scaff) + 8))
#define GET_CONTIG_ID(contig) ((contig)[6] == '_' ? atoi((contig) + 7) : atoi((contig) + 6))

#define GET_READ(s, i) ((s) + ((_cfg->max_read_len + 1) * (i)))

const char *reject_reasons[] = {
    "FORMAT", "MINLEN", "UNINFORMATIVE", "5-TRUNCATED", "3-TRUNCATED",
    "SINGLETON", "NO_ASSOCIATED_GAP", "TRUNC5_FAIL", "TRUNC3_FAIL", 
    "NO_SCAFFOLD_INFO", "POTENTIAL_INNIE", "POTENTIAL_SHORTY", "FILE HEADER"};

//#define FTIME_TO_LINE(t) tprintf_flush("%s, line %d: time %.4f\n", __func__, __LINE__, ELAPSED_TIME(t))

typedef struct {
    //char name[MAX_READ_NAME_LEN];
    int q_start, q_stop, q_length;
    int s_start, s_stop, s_length;
    int contig_id;
    char strand;
    //int score;
    //int e_value;
    //int identities;
    //int align_length;
    char start_status, end_status;
    int dirn;
} alignment_t;

typedef struct {
    char read_name[MAX_READ_NAME_LEN];
    int dirn;
    int gap_id;
    char strand;
    int contig_id;
} read_gap_pair_t;

typedef struct {
    char name[MAX_READ_NAME_LEN];
    long aln;
} read_aln_t;

typedef struct {
    char name[MAX_READ_NAME_LEN];
    int gap_id;
    char orient;
} read_orient_t;

typedef struct {
    read_orient_t *reads;
    long num;
    long size;
} read_list_t;

typedef struct {
    // if we store an index into a separate array of alignments instead of the actual 
    // alignments, we can probably reduce space usage by 50%
    // however, the code is more complicated, requiring two different structures instead
    // of one. It's only worth doing if we are running out of memory, or getting too 
    // high a loss of excessive alns.
    alignment_t alns[MAX_PAIR_ALNS];
    char name[MAX_READ_NAME_LEN];
    int nalns;
} pair_alns_t;

typedef struct {
    shared [] pair_alns_t *pairs;
    long npairs;
    long next_pair;
} pair_alns_array_t;

DEFINE_HTABLE_TYPE(read_gap_pair);
DEFINE_HTABLE_TYPE(read_aln);

static config_t *_cfg;

static int _ngaps_block = 0;
static int _max_ngaps = 0;
static int _tot_ngaps = 0;
static shared contig_t *_contigs = NULL;
static int _max_contigs = 0;
static shared scaffold_t *_scaffolds;
static int _max_scaffolds = 0;
static shared gap_set_t *_gaps = NULL;
static int _max_contigs_per_scaffold = 0;
static long _nreads = 0;

static inline gap_t *get_mygap(int gi)
{
    int thread = gi / _ngaps_block;
    if (gi < _gaps[MYTHREAD].gap_start || gi >= _gaps[MYTHREAD].gap_end || thread != MYTHREAD)
        DIE("gap %d is not in our gaps\n", gi);
    int gapi = gi % _ngaps_block;
    if (gapi < 0 || gapi >= _gaps[thread].ngaps)
        DIE("Gapi out of range in get_gapi: %d ! in [%d, %d)\n", gapi, 0, _gaps[thread].ngaps);
    return (gap_t*)&(_gaps[thread].gaps[gapi]);
}

// this is also called from merauder, hence not inline
shared [] gap_t *get_gap(int gi, int line, const char *fname, int abort_on_fail)
{
    int thread = gi / _ngaps_block;
    if (thread < 0 || thread >= THREADS) 
        DIE("At %s:%d: Thread out of range in get_gapi: %d\n", fname ? fname : __FILE__, line, thread);
    int gapi = gi % _ngaps_block;
    if (gapi < 0 || gapi >= _gaps[thread].ngaps) {
        if (abort_on_fail == ABORT_ON_FAIL)
            DIE("At %s:%d, Gapi out of range in get_gapi: %d ! in [%d, %d)\n", 
                fname ? fname : __FILE__, line, gapi, 0, _gaps[thread].ngaps);
        else
            return NULL;
    }
    return &_gaps[thread].gaps[gapi];
}

char *get_gapi_info(int gi, char *buf)
{
    shared [] gap_t *gap = get_gap(gi, __LINE__, NULL, ABORT_ON_FAIL);
    sprintf(buf, "Scaffold%d:Contig%d.%d<-[%d +/- %.0f]->Contig%d.%d", 
            gap->scaff_id, gap->contig1_id, gap->contig_ext1,
            gap->size, gap->uncertainty, gap->contig2_id,
            gap->contig_ext2);
    return buf;
}

shared contig_t *get_contig(int id)
{
    if (id > _max_contigs || id < 0) {
        //WARN("Requested contig out of range %d > %d\n", id, _max_contigs);
        return NULL;
    }
    assert(_contigs && IS_VALID_UPC_PTR(_contigs));
    if (_contigs[id].id == -1) {
        //WARN("No Contig%d exists\n", id);
        return NULL;
    }
    return &_contigs[id];
}

shared scaffold_t *get_scaffold(int id)
{
    if (id > _max_scaffolds || id < 0)
        DIE("Requested scaffold out of range %d >= %d\n", id, _max_scaffolds);
    if (_scaffolds[id].id == -1)
        DIE("No scaffold%d exists\n", id);
    return &_scaffolds[id];
}

shared contig_t *get_contigs(void)
{ 
    assert(_contigs && IS_VALID_UPC_PTR(_contigs));
    return _contigs;
}

shared scaffold_t *get_scaffolds(void) 
{
    assert(_scaffolds && IS_VALID_UPC_PTR(_scaffolds));
    return _scaffolds;
}

int get_max_contigs(void)
{
    return _max_contigs;
}

int get_max_scaffolds(void)
{
    return _max_scaffolds;
}

void init_globals(config_t *cfg, shared gap_set_t *gaps) 
{
    _cfg = cfg;
    _gaps = gaps;
#ifdef DBG_GET_ALIGNMENT
    serial_printf("DBG_GET_ALIGNMENT\n");
#endif
#ifdef DBG_PAIR_ALNS
    serial_printf("DBG_PAIR_ALNS\n");
#endif
#ifdef DBG_PUT_READ_ORIENT
    serial_printf("DBG_PUT_READ_ORIENT\n");
#endif
#ifdef DBG_READ_DATA
    serial_printf("DBG_READ_DATA\n");
#endif
#ifdef DBG_CONTIGS
    serial_printf("DBG_CONTIGS\n");
#endif    
#ifdef DBG_SEQUENCES
    serial_printf("DBG_SEQUENCES\n");
#endif
#ifdef DBG_SCAFFOLDS
    serial_printf("DBG_SCAFFOLDS\n");
#endif
}

char *get_report_line(gap_t *gap, Buffer buf)
{
    printfBuffer(buf, "REPORT:\tScaffold%d\tContig%d.%d\t%s\tContig%d.%d\t%s\t%d\t%.0f",
             gap->scaff_id, gap->contig1_id, gap->contig_ext1, 
             gap->primer1, gap->contig2_id, gap->contig_ext2,
             gap->primer2, gap->size, gap->uncertainty);
    return buf->buf;
}

char *get_gap_info(gap_t *gap, char *buf, int nchar)
{
    snprintf(buf, nchar, "Scaffold%d:Contig%d.%d<-[%d +/- %.0f]->Contig%d.%d", 
             gap->scaff_id, gap->contig1_id, gap->contig_ext1,
             gap->size, gap->uncertainty, gap->contig2_id,
             gap->contig_ext2);
    return buf;
}


static char *get_aln_info(alignment_t *aln, char *name, char *buf)
{
    sprintf(buf, "'%s'/%d\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%c\t%c\t%c", 
            name, aln->dirn, aln->q_start, aln->q_stop, aln->q_length, 
            aln->contig_id, aln->s_start, aln->s_stop, aln->s_length, aln->strand,
            aln->start_status, aln->end_status);
    return buf;
}


void fprint_gap(FILE *f, gap_t *gap)
{
    fprintf(f, "Scaffold%d\tContig%d.%d\t%s\tContig%d.%d\t%s\t%d\t%.0f", 
            gap->scaff_id, gap->contig1_id, gap->contig_ext1, 
            gap->primer1, gap->contig2_id, gap->contig_ext2, 
            gap->primer2, gap->size, gap->uncertainty);
    for (int j = 0; j < gap->nreads; j++) 
        fprintf(f, "\t%s:%s", (char*)GET_READ(gap->reads_nts, j), (char*)GET_READ(gap->reads_quals, j));
    fprintf(f, "\n");
}

#ifdef SHOW_MEMORY
void print_memory_used(void) 
{
    memory_usage_t mem;
    get_mem_usage(&mem);
	tprintf_flush("Memory used (GB): VmSize %.2f, VmRSS %.2f, share %.2f\n",
                  (double)mem.vm_size / ONE_GB, (double)mem.vm_rss / ONE_GB, 
                  (double)mem.share / ONE_GB);
}
#else
#define print_memory_used()
#endif 

#ifndef HIPMER_MAX_DEPTH
#define HIPMER_MAX_DEPTH 20000
#endif
#define MAX_DEPTH HIPMER_MAX_DEPTH 
#define DEPTH_HIST(thread, depth) depth_hist[(thread) + THREADS * (depth)]

static shared int _peak_depth = 0;

static void calc_scaffold_depths(void)
{
    shared double *depth_hist = (shared double *) upc_all_alloc(MAX_DEPTH * THREADS, sizeof(double));
    if (depth_hist == NULL) DIE("Could not allocate %lld bytes over %d threads\n", (long long) THREADS * sizeof(double) * MAX_DEPTH, THREADS);
    for (int depth = 0; depth < MAX_DEPTH; depth++) {
        assert( upc_threadof(&DEPTH_HIST(MYTHREAD,depth)) == MYTHREAD );
        assert( (&DEPTH_HIST(MYTHREAD, depth) - depth_hist) < MAX_DEPTH * THREADS);
        DEPTH_HIST(MYTHREAD, depth) = 0;
    }

    upc_forall (int scaffidx = 0; scaffidx < _max_scaffolds; scaffidx++; scaffidx) {
        if (_scaffolds[scaffidx].id == -1) {
            continue;
        }
        double mean_depth = _scaffolds[scaffidx].depth / _scaffolds[scaffidx].weight;
        _scaffolds[scaffidx].depth = mean_depth;
        if ((int)mean_depth < 0) 
            DIE("mean_depth %f < 0\n", mean_depth);
        int md = bankers_round(mean_depth);
        if (md < 0 || md >= MAX_DEPTH) {
            WARN("Depth hist is too small for scaffold %d: %d < %d. scaffolds.depth=%f scaffolds.weight=%f\n", 
                 _scaffolds[scaffidx].id, MAX_DEPTH, md, _scaffolds[scaffidx].depth,  _scaffolds[scaffidx].weight);
        } else {
            assert( (&DEPTH_HIST(MYTHREAD, md) - depth_hist) < MAX_DEPTH * THREADS);
            DEPTH_HIST(MYTHREAD, md) += _scaffolds[scaffidx].weight;
        }
    }
    UPC_TIMED_BARRIER;
    upc_forall (int depth = 0; depth < MAX_DEPTH; depth++; depth) {
        double depthCount = 0;
        for (int thread = 0; thread < THREADS; thread++) {
            depthCount += DEPTH_HIST(thread, depth);
        }
        DEPTH_HIST(0,depth) = depthCount;
    }
    UPC_TIMED_BARRIER;
    if (!MYTHREAD) {
	double max_weight = 0;
	for (int depth = 0; depth < MAX_DEPTH; depth++) {
            if (DEPTH_HIST(0,depth) > max_weight) {
                _peak_depth = depth;
                assert( upc_threadof( &(DEPTH_HIST(0,depth)) ) == MYTHREAD );
                max_weight = DEPTH_HIST(0,depth);
            }
        }
    }
    UPC_TIMED_BARRIER;
    if (!MYTHREAD)
        upc_free(depth_hist);

    tprintf_flush("Modal scaffold depth %d\n", _peak_depth);
    upc_forall (int i = 0; i < _max_scaffolds; i++; i) {
        if (_scaffolds[i].id == -1) {
            continue;
        }
        assert(_scaffolds && _scaffolds[i].depth >= 0.0);
        _scaffolds[i].depth /= _peak_depth;
#ifdef DBG_SCAFFOLDS
        for (int j = 0; j < _scaffolds[i].ncontigs; j++) {
            dbg("DSCAF %d\t%d\t%.0f\t%.0f\t%.0f\n", 
                _scaffolds[i].id, _scaffolds[i].contig_ids[j],
                _scaffolds[i].length, _scaffolds[i].depth, _scaffolds[i].weight);
        }
#endif
    }

    UPC_TIMED_BARRIER;
}

static void get_srf_counts(FILE *f, char *fname)
{
    // read file to get gap count
    _gaps[MYTHREAD].ngaps = 0;
	int max_scaff_id = 0;
	int max_contig_id = 0;
    int max_contigs_per_scaffold = 0;
    int contigs_per_scaffold = 1;
    int slen = strlen("Scaffold");
    int clen = strlen("CONTIG");
    const int buf_size = 256;
    char buf[256];
    int prev_scaff_id = -1;
    while (fgets(buf, buf_size, f)) {
        CHECK_LINE_LEN(buf, buf_size, fname);
        if (strstr(buf, "GAP")) {
			_gaps[MYTHREAD].ngaps++;
		} else {
			char *end_scaff = strchr(buf, '\t');
			if (!end_scaff)
				DIE("couldn't get tab\n");
			end_scaff[0] = 0;
			int scaff_id = atoi(buf + slen);
			if (scaff_id > max_scaff_id)
				max_scaff_id = scaff_id;
			if (prev_scaff_id == scaff_id) {
				contigs_per_scaffold++;
			} else {
				if (contigs_per_scaffold > max_contigs_per_scaffold) {
					max_contigs_per_scaffold = contigs_per_scaffold;
				}
				contigs_per_scaffold = 1;
				prev_scaff_id = scaff_id;
			}
			char *contig = strstr(end_scaff + 1, "Contig");
			if (!contig)
				DIE("Couldn't find the Contig in '%s'\n", buf);
			char *end_contig = strchr(contig, '\t');
			if (!end_contig)
				DIE("couldn't get tab\n");
			end_contig[0] = 0;
			int contig_id = atoi(contig + clen);
			if (contig_id > max_contig_id)
				max_contig_id = contig_id;
		}
    }
    if (contigs_per_scaffold > max_contigs_per_scaffold) {
		max_contigs_per_scaffold = contigs_per_scaffold;
    }
    // reduce to max
    _ngaps_block = reduce_int(_gaps[MYTHREAD].ngaps, UPC_MAX);
    _max_ngaps = _ngaps_block * THREADS;
    _tot_ngaps = reduce_int(_gaps[MYTHREAD].ngaps, UPC_ADD);

    tprintf_flush("Max gaps in block %d, max gaps total %d\n", _ngaps_block, _tot_ngaps);
    tprintf_flush("Max contig id %d, max scaff id %d\n", max_contig_id, max_scaff_id);
    tprintf_flush("Max contigs per scaffold %d\n", max_contigs_per_scaffold);

    _max_contigs = reduce_int(max_contig_id + 1, UPC_MAX);
    _max_scaffolds = reduce_int(max_scaff_id + 1, UPC_MAX);
    _max_contigs_per_scaffold = reduce_int(max_contigs_per_scaffold + 1, UPC_MAX);

    tprintf_flush("Shared allocation: %d contigs and %d scaffolds\n", _max_contigs, _max_scaffolds);
    rewind(f);
}

// Reads in scaffolds and computes weight and depth statistics
// Only reads those scaffolds relevant to this thread (set of gaps)
// Format of file:
// scaffold_id CONTIGi contig_id start end depth
// scaffold_id GAPi gap_size gap_uncertainty
void process_scaffold_report_file(void)
{
    START_TIMER(T_PROCESS_SRF);
    const int buf_size = 256;
    char *buf = malloc_chk(buf_size);

    char fname[MAX_FILE_PATH];
    sprintf(fname, "%s/%s_%d", _cfg->base_dir, _cfg->sr_file, MYTHREAD);
    get_rank_path(fname, MYTHREAD);
    FILE *f = fopen_chk(fname, "r");
    START_TIMER(T_SRF_COUNTS);
    get_srf_counts(f, fname);
    stop_timer(T_SRF_COUNTS);

    // allocate at least the size of the block, so we can balance later
    int ngaps = _gaps[MYTHREAD].ngaps;
    if (ngaps < _ngaps_block)
        ngaps = _ngaps_block;
    _gaps[MYTHREAD].gaps = upc_alloc(ngaps * sizeof(gap_t));
    _gaps[MYTHREAD].gap_start = MYTHREAD * _ngaps_block;
    _gaps[MYTHREAD].gap_end = _gaps[MYTHREAD].gap_start + _gaps[MYTHREAD].ngaps;

    for (int i = 0; i < _gaps[MYTHREAD].ngaps; i++) {
        gap_t *gap = (gap_t*)&(_gaps[MYTHREAD].gaps[i]);
        gap->nreads = 0;
        gap->max_reads = 0;
    }
    
    _contigs = upc_all_alloc(_max_contigs, sizeof(contig_t));
    if (_contigs == NULL) DIE("Could not allocate %d contigs\n", _max_contigs);
    upc_forall (int i = 0; i < _max_contigs; i++; i) {
        _contigs[i].id = -1;
        _contigs[i].seq = NULL;
        _contigs[i].seq_len = 0;
        _contigs[i].gapi_3 = -1;  // undef
        _contigs[i].gapi_5 = -1;
    }

    _scaffolds = upc_all_alloc(_max_scaffolds, sizeof(scaffold_t));
    if (_scaffolds == NULL) DIE("Could not allocate %d scaffolds\n", _max_scaffolds);
    upc_forall (int i = 0; i < _max_scaffolds; i++; i) {
        _scaffolds[i].id = -1;
        _scaffolds[i].ncontigs = 0;
        _scaffolds[i].depth = 0;
        _scaffolds[i].weight = 0;
        _scaffolds[i].contig_ids = upc_alloc(_max_contigs_per_scaffold * sizeof(int));
    }
    UPC_TIMED_BARRIER;

    const int max_toks = 6;
    char **tokens = malloc_chk(max_toks * sizeof(char *));
    int tot_len;
    shared scaffold_t *scaffold = NULL;
    shared contig_t *contig = NULL;
    int contig_id = -1;
    int scaff_start = 0;
    int scaff_end = 0;
    int open_gap = 0;
    int scaff_depths_avail = 0;
    int nscaffolds = 0;
    int ncontigs = 0;
    int gapi = 0;
    serial_printf("Processing scaffold_report %s... ", fname);
    tprintf_flush("Processing scaffold_report %s... ", fname);
    long line = -1;
    START_TIMER(T_SRF_READS);
    while (fgets(buf, buf_size, f)) {
        line++;
        int ntoks = get_tokens(buf, tokens, &tot_len, max_toks);
        if (tokens[1][0] == 'C') {
            int scaff_id = GET_SCAFF_ID(tokens[0]);
            if (scaff_id > _max_scaffolds || scaff_id < 0)
                DIE("scaffold id %d is too high > %d\n", scaff_id, _max_scaffolds);
            assert(_scaffolds && IS_VALID_UPC_PTR(_scaffolds));
            if (_scaffolds[scaff_id].id == -1) {
                _scaffolds[scaff_id].id = scaff_id;
                nscaffolds++;
            }
            scaffold = &_scaffolds[scaff_id];

            contig_id = GET_CONTIG_ID(tokens[2] + 1);
            if (contig_id >= _max_contigs || contig_id < 0)
                DIE("contig id %d is too high > %d\n", contig_id, _max_contigs);
            if (_contigs[contig_id].id != -1) 
                DIE("Current version only supports unique contig->scaffold mappings: "
					 "Contig%d, Scaffold%d, previous Contig%d\n", 
					 contig_id, scaff_id, _contigs[contig_id].id);

            assert(_contigs && IS_VALID_UPC_PTR(_contigs));
            contig = &_contigs[contig_id];
            contig->id = contig_id;
            contig->strand = tokens[2][0];
            contig->s_start = atoi(tokens[3]);
            contig->s_end = atoi(tokens[4]);
            ncontigs++;

            int contig_i = bupc_atomicI_fetchadd_relaxed(&(scaffold->ncontigs), 1);
            if (contig_i >= _max_contigs_per_scaffold || contig_i < 0) 
                DIE("Scaffold %d - Contig id out of range: %d >= %d\n", 
                     scaff_id, contig_i, _max_contigs_per_scaffold);

            scaffold->contig_ids[contig_i] = contig_id;

            scaff_start = contig->s_start;
            scaff_end = contig->s_end;
            scaffold->length = atoi(tokens[4]);
            if (ntoks > 5) {
                scaff_depths_avail = 1;
                double depth = atof(tokens[5]);
                double weight = scaff_end - scaff_start + 1;
                scaffold->depth += (depth * weight);
                scaffold->weight += weight;
            } else if (_cfg->exclude_repeats > 0) {
                tprintf("Warning: Depth information not available. "
                        "Repeat exclusion turned off.\n");
                _cfg->exclude_repeats = 0;
            }
            if (open_gap) {
                int gi = (gapi - 1) + _gaps[MYTHREAD].gap_start;
                gap_t *gap = get_mygap(gi);
                if (gi  != gap->id)
                    DIE("Incorrect gap id: expected %d, got %d\n", gi, gap->id);
                gap->contig2_id = contig->id;
                if (contig->strand == '+') 
                    gap->contig_ext2 = 5;
                else 
                    gap->contig_ext2 = 3;
                if (contig->strand == '+') 
                    contig->gapi_3 = gi;
                else 
                    contig->gapi_5 = gi;
                open_gap = 0;
            }
        } else { // a gap
            if (tokens[1][0] != 'G')
                DIE("Expected a GAP at line %ld in %s\n", line, fname);
            int gi = gapi + _gaps[MYTHREAD].gap_start;
            if (!contig)
                DIE("couldn't get contig for gap %d\n", gi);
            gap_t *gap = get_mygap(gi);
            gap->id = gi;
            gap->size = atoi(tokens[2]);;
            gap->uncertainty = atoi(tokens[3]);
            if (!scaffold)
                DIE("couldn't get scaffold for gap %d\n", gi);
            gap->scaff_id = scaffold->id;
            gap->start = scaff_end;
            gap->end = gap->start + gap->size + 1;
            if (contig->strand == '+') 
                gap->contig_ext1 = 3;
            else 
                gap->contig_ext1 = 5;
            gap->contig1_id = contig->id;
            gap->nreads = 0;
            if (contig->strand == '+') 
                contig->gapi_5 = gi;
            else 
                contig->gapi_3 = gi;
            open_gap = 1;
            gapi++;
        }
    }
    stop_timer(T_SRF_READS);
    fclose(f);
    free(buf);
    free(tokens);
    if (gapi != _gaps[MYTHREAD].ngaps)
        DIE("Wrong gap count second time around %d != %d\n", gapi, _gaps[MYTHREAD].ngaps);
    tprintf_flush("\n... %d gaps found, %d scaffolds found, %d contigs found\n", 
                  _gaps[MYTHREAD].ngaps, nscaffolds, ncontigs);
    UPC_TIMED_BARRIER;
#ifdef DBG_CONTIGS
    char gbuf[256];
    assert(_contigs && IS_VALID_UPC_PTR(_contigs));
    for (int i = 0; i < _max_contigs; i++) {
        if (_contigs[i].id == -1) {
            continue;
        }
        dbg("DCONT %cContig%d\t%d\t%d\n", _contigs[i].strand, _contigs[i].id, 
            _contigs[i].s_start, _contigs[i].s_end);
        /*
		  dbg("DCONT %d\t%c\t", _contigs[i].id, _contigs[i].strand);
		  if (_contigs[i].gapi_3 != -1)
		  dbg("%s\t", get_gapi_info(_contigs[i].gapi_3, gbuf));
		  else 
		  dbg("-1\t");
		  if (_contigs[i].gapi_5 != -1)
		  dbg("%s\n", get_gapi_info(_contigs[i].gapi_5, gbuf));
		  else 
		  dbg("-1\n");
        */
    }
#endif
    START_TIMER(T_CALC_SCAFF_DEPTH);
    scaff_depths_avail = reduce_int(scaff_depths_avail, UPC_MAX); 
    if (scaff_depths_avail) {
        calc_scaffold_depths();
    }
    stop_timer(T_CALC_SCAFF_DEPTH);
    stop_timer(T_PROCESS_SRF);
    serial_printf("done in %.3f s\n", get_elapsed_time(T_PROCESS_SRF));
    serial_printf("Modal scaffold depth %d\n", _peak_depth);
    print_memory_used();
}

// get list of contigs and add to hash table of contigs
// file format:
// >contig_id
// contig (split over multiple lines)
//
// Each thread reads the same whole contigs file. It only gets the data it needs
// for the gaps it processes. This makes the read quicker and uses less memory.
// An attempt to do this as a shared data structure, loading the whole file into
// memory is actually a lot slower. So we'll keep it this simpler way.
void process_contigs_file(void)
{
    START_TIMER(T_PROCESS_CONTIGS);

    assert(_contigs && IS_VALID_UPC_PTR(_contigs));
    long num_entries = 0;
    shared contig_t *contig;
    const int buf_size = 100;
    char fname[MAX_FILE_PATH];
    sprintf(fname, "%s/%s_%d.fasta", _cfg->base_dir, _cfg->contigs_file, MYTHREAD);
    get_rank_path(fname, MYTHREAD);
    serial_printf("Processing %s... ", fname);
    tprintf_flush("Processing %s... ", fname);
    gzFile f = gzopen(fname, "r");
    if (!f) { DIE("Could not open %s\n", fname); }
    kseq_t *ks = kseq_init(f);
    while ( kseq_read(ks)>=0 ) {
        int contig_id = GET_CONTIG_ID(ks->name.s);
        if (contig_id >= _max_contigs || contig_id < 0) {
            tprintf_flush("Skipping unused contig: contig_id %d > max_contigs %d\n", 
                          contig_id, _max_contigs);
            continue;
        }
        assert(contig_id >= 0);
        assert(_contigs && IS_VALID_UPC_PTR(_contigs));
        if (_contigs[contig_id].id == -1) {
            continue;
        }
        contig = &_contigs[contig_id];
        int seq_len = ks->seq.l;
        contig->seq = upc_alloc(seq_len + 1);
        if (!contig->seq) DIE("Could not allocate %d for new seq\n", seq_len);
        upc_memput(contig->seq, ks->seq.s, seq_len + 1);
        contig->seq_len = seq_len;
        num_entries++;
    }
    kseq_destroy(ks);
    stop_timer(T_PROCESS_CONTIGS);
    serial_printf("done in %.3f s\n", get_elapsed_time(T_PROCESS_CONTIGS));
    tprintf_flush("\n... %d contigs added\n", num_entries);
    print_memory_used();
    UPC_TIMED_BARRIER;
    //double t_all = get_timer_all_max(T_PROCESS_CONTIGS);
    //serial_printf("All done in %.3f s\n", t_all);
#ifdef DBG_SEQUENCES
    if (MYTHREAD == 0) {
        assert(_contigs && IS_VALID_UPC_PTR(_contigs));
        for (int i = 0; i < _max_contigs; i++) {
            if (_contigs[i].id != -1) {
                UPC_MEMGET_STR(seq, _contigs[i].seq, _contigs[i].seq_len + 1);
                if  (_contigs[i].id != i)
                    DIE("wrong contig id %d != %d\n", _contigs[i].id, i);
                dbg("DSEQ %d %d %s\n", _contigs[i].id, _contigs[i].seq_len, seq);
            }
        }
    }
#endif
}

static inline char classify_alignment(int unaligned, int projected, int s_length,
                                      int wiggle_room, int truncate, int truncate_val) 
{
    int projected_off = 0;
    if (projected < 1)
        projected_off = 1 - projected;
    else if (projected > s_length)
        projected_off = projected -  s_length;
    int missing_bases = unaligned - projected_off;
    // classify this alignment
    if (unaligned == 0) 
        return 'F'; // FUL
    if (projected_off > 0 && missing_bases < wiggle_room) 
        return 'G'; // GAP
    if (unaligned < wiggle_room || truncate == truncate_val) 
        return 'I'; // INC
    return 0; // rejected
}


#define MAX_ALIGNMENT_TOKENS 14
static int get_alignment(char *buf, alignment_t *aln, char *pair_name)
{
    assert(buf != NULL);
    assert(aln != NULL);
    START_TIMER(T_GET_ALN);
    // 14 fields, tab separated
    const int max_toks = MAX_ALIGNMENT_TOKENS;
    char *tokens[MAX_ALIGNMENT_TOKENS];
    int tot_len;
    int ntoks = get_tokens(buf, tokens, &tot_len, max_toks);
    if (ntoks != max_toks)
        return REJECT_FORMAT;
    assert(strlen(tokens[0]) > 11);
    if (tokens[0][11] == 'F')
        return REJECT_NO_ASSOCIATED_GAP;
    //aln->identities = atoi(tokens[12]);
    //if (aln->identities < _cfg->mer_size) 
    //   return REJECT_MINLEN;
    aln->contig_id = GET_CONTIG_ID(tokens[5]);
    shared contig_t *contig = get_contig(aln->contig_id);
    if (!contig)
        return REJECT_NO_SCAFFOLD_INFO;
    //get_contig_chk(contig, aln->contig_id);
    // If there's no gap associated with the contig, skip the alignment
    if (contig->id == -1)
        return REJECT_NO_ASSOCIATED_GAP;
    if (contig->gapi_3 == -1 && contig->gapi_5 == -1)
        return REJECT_NO_ASSOCIATED_GAP;
    char *name;
    if (!get_fq_name_dirn(tokens[1], &name, &aln->dirn))
        DIE("Error: couldn't extract pairing info from the read name:  %s\n", tokens[1]);
    if (!pair_name[0])
        strcpy(pair_name, name);

    aln->q_start = atoi(tokens[2]);
    aln->q_stop = atoi(tokens[3]);
    aln->q_length = atoi(tokens[4]);
    
    aln->s_start = atoi(tokens[6]);
    aln->s_stop = atoi(tokens[7]);
    aln->s_length = atoi(tokens[8]);
    char *strand = tokens[9];
    if (strcmp(strand, "Plus") != 0 && strcmp(strand, "Minus") != 0)
        DIE("Invalid entry in meraligner file:\n%s\n", buf);
    aln->strand = (strand[0] == 'P' ? '+' : '-');
    //aln->score = atoi(tokens[10]);
    //aln->e_value = atoi(tokens[11]);

    //aln->align_length = atoi(tokens[13]);
    aln->start_status = 0;
    aln->end_status = 0;
    stop_timer(T_GET_ALN);
#ifdef DBG_GET_ALIGNMENT
    char alnbuf[200];
    dbg("DGET %s\n", get_aln_info(aln, name, alnbuf));
#endif
    return GOOD_ALIGNMENT;
}

static void add_read_gap_pair(htable_t pairs, int gap_id, alignment_t *aln, char *name)
{
    char buf[100];
    sprintf(buf, "%s/%d:%d", name, aln->dirn, gap_id);
    read_gap_pair_t *pair = htable_get_read_gap_pair(pairs, buf);
    if (!pair) {
        pair = malloc_chk(sizeof(read_gap_pair_t));
        strcpy(pair->read_name, name);
        pair->gap_id = gap_id;
        CHECK_ERR(htable_put_read_gap_pair(pairs, strdup(buf), pair, NO_CHECK_DUPS));
    }
    pair->dirn = aln->dirn;
    pair->strand = aln->strand;
    pair->contig_id = aln->contig_id;
}

// to compute alignment, sum up 3 constants, e.g. LEFT+PLUS+FIVE means the read
// is aligned on the Plus strand to the contig to the Left of the gap whose 5'
// end is adjacent to the gap
static const int LEFT = 0;
static const int PLUS = 0;
static const int FIVE = 0;
static const int THREE = 1;
static const int MINUS = 2;
static const int RIGHT = 4;
static const char combos[] = {'-','+','+', '-', '+', '-', '-', '+'};

static char get_orient(int contig_id, char strand, int gi)
{
    if (gi == -1)
        DIE("gi is -1\n");
    int combo = 0;
    shared [] gap_t *gap = get_gap(gi, __LINE__, NULL, ABORT_ON_FAIL);
    int cid_1 = gap->contig1_id;
    int cid_2 = gap->contig2_id;
    if (cid_1 == contig_id) 
        combo = LEFT + (gap->contig_ext1 == 3 ? THREE : FIVE);
    else if (cid_2 == contig_id)
        combo = RIGHT + (gap->contig_ext2 == 3 ? THREE : FIVE);
    else {
        DIE("Neither Contig%d nor Contig%d match Contig%d for gap %d\n", 
             cid_1, cid_2, contig_id, gi);
    }
    combo += (strand == '+' ? PLUS : MINUS);
    return combos[combo];
}

static int _num_reallocs = 0;

static void add_to_read_list(read_list_t *read_list, char *read_name, int gap_id, char orient)
{
    read_orient_t *read_orient = &read_list->reads[read_list->num];

    strcpy(read_orient->name, read_name);
    read_orient->gap_id = gap_id;
    read_orient->orient = orient;
    read_list->num++;
    if (read_list->num == read_list->size) {
        read_list->size *= 1.5;
        read_list->reads = realloc_chk(read_list->reads, read_list->size * sizeof(read_orient_t));
        _num_reallocs++;
    }
    // increment the number of reads for this gap
    bupc_atomicI_fetchadd_relaxed(&(get_gap(gap_id, __LINE__, NULL, ABORT_ON_FAIL)->max_reads), 1);
}

static long direct_alignments(pair_alns_t *pair_alns, read_list_t *read_list)
{
    START_TIMER(T_DIRECT_ALNS);
    const long max_pairs = pair_alns->nalns * 2;
    htable_t pairs = create_htable(max_pairs * 4, "pairs");
    long nplaced = 0;
    for (long i = 0; i < pair_alns->nalns; i++) {
        alignment_t *aln = &(pair_alns->alns[i]);
        shared contig_t *contig;
        get_contig_chk(contig, aln->contig_id);
        // If the read aligns into a gap, add a read<->gap pair
        if (aln->start_status == 'G') {
            int gap_id = aln->strand == '+' ? contig->gapi_3 : contig->gapi_5;
            if (gap_id != -1)
                add_read_gap_pair(pairs, gap_id, aln, pair_alns->name);
        }
        if (aln->end_status == 'G') {
            int gap_id = aln->strand == '+' ? contig->gapi_5 : contig->gapi_3;
            if (gap_id != -1)
                add_read_gap_pair(pairs, gap_id, aln, pair_alns->name);
        }
    }
    char read_name[MAX_READ_NAME_LEN];
    htable_iter_t iter = htable_get_iter(pairs);
    read_gap_pair_t *pair;
    // For reads aligning directly into a gap, record the orientation of the
    // read wrt to the gap
    while ((pair = htable_get_next_read_gap_pair(pairs, iter)) != NULL) {
        sprintf(read_name, "%s/%d", pair->read_name, pair->dirn);
        char orient = get_orient(pair->contig_id, pair->strand, pair->gap_id);
        // Store the read_name <-> (gap_id, orient) mapping. We then later use
        // this read_name as a lookup for getting nts and qual info from the
        // fastq file
        add_to_read_list(read_list, read_name, pair->gap_id, orient);
        nplaced++;
        free(pair);
    }
    free(iter);
    destroy_htable(pairs, FREE_KEYS);
    stop_timer(T_DIRECT_ALNS);
    return nplaced;
}

static char *get_mate_name(char *name)
{
    char *mate = strdup(name);
    int mi = strlen(mate) - 1;
    mate[mi] = (name[mi] == '1' ? '2' : '1');
    return mate;
}

static long projected_alignments(pair_alns_t *pair_alns, read_list_t *read_list, int insert_size, 
                                 int insert_sigma)
{
    START_TIMER(T_PROJECTED_ALNS);
    const long max_pairs = pair_alns->nalns * 2;
    htable_t reads_to_alns = create_htable(max_pairs * 4, "reads_to_alns");
    int max_aligned_scaff_len = 0;
    long nprojected = 0;
    for (long i = 0; i < pair_alns->nalns; i++) {
        alignment_t *aln = &(pair_alns->alns[i]);
        shared contig_t *contig;
        get_contig_chk(contig, aln->contig_id);
        // one alignment is recorded for each read to be used to project the
        // read's mate if needed the alignment to the largest scaffold is
        // retained 
        // each read is assigned to a single contig/scaffold for this purpose
        // mate pairs may only be projected into gap(s) on a single scaffold by
        // this construction (although direct read alignments may place reads in
        // gaps on multiple scaffolds)
        int gap_id = aln->strand == '+' ? contig->gapi_5 : contig->gapi_3;
        if (gap_id == -1)
            continue;
        shared scaffold_t *scaffold = get_scaffold(get_gap(gap_id, __LINE__, NULL, ABORT_ON_FAIL)->scaff_id);
        int scaff_len = scaffold->length;
        char buf[MAX_READ_NAME_LEN];
        sprintf(buf, "%s/%d", pair_alns->name, aln->dirn);
        read_aln_t *read_aln = htable_get_read_aln(reads_to_alns, buf);
        if (read_aln) {
            if (scaff_len > max_aligned_scaff_len) {
                read_aln->aln = i;
                max_aligned_scaff_len = scaff_len;
            }
        } else {
            read_aln = malloc_chk(sizeof(read_aln_t));
            strcpy(read_aln->name, buf);
            read_aln->aln = i;
            CHECK_ERR(htable_put_read_aln(reads_to_alns, read_aln->name, read_aln, NO_CHECK_DUPS));
            max_aligned_scaff_len = 0;
        }
    }
    const int project_z = 2;
    int max_project = insert_size + project_z * insert_sigma;
    int min_project = insert_size - project_z * insert_sigma;
    read_aln_t *read_aln;
    htable_iter_t iter = htable_get_iter(reads_to_alns);
    while ((read_aln = htable_get_next_read_aln(reads_to_alns, iter)) != NULL) {
        char *mate_name = get_mate_name(read_aln->name);
        if (!mate_name)
            DIE("Could not determine mate name for read %s\n", read_aln->name);
        // Attempt to project the unplaced mate into a gap
        // assumes unplaced mate is same length as the placed one
        if (htable_get_read_aln(reads_to_alns, mate_name)) {
            free(mate_name);
            continue;
        }
        // only project if the mate has no alignment
        alignment_t *aln = &(pair_alns->alns[read_aln->aln]);
        int unaligned_start = aln->q_start - 1;
        int far_projection = max_project - unaligned_start;
        int near_projection = min_project - aln->q_length - unaligned_start;
        int start_gap = -1;
        shared contig_t *contig;
        get_contig_chk(contig, aln->contig_id);
        if (aln->strand == '+') {
            far_projection = aln->s_start + far_projection;
            near_projection = aln->s_start + near_projection;
            start_gap = contig->gapi_5;
        } else {
            far_projection = aln->s_stop - far_projection;
            near_projection = aln->s_stop - near_projection;
            start_gap = contig->gapi_3;
        }
        if (start_gap == -1) {
            free(mate_name);
            continue;
        }
        int test_gapi = start_gap;
        shared [] gap_t *sgap = get_gap(start_gap, __LINE__, NULL, ABORT_ON_FAIL);
        int scaff_id = sgap->scaff_id;
        int cid_1 = sgap->contig1_id;
        int cid_2 = sgap->contig2_id;
        // find anchor-contig origin position in scaffold coord, assign
        // direction of scan via iterator and transform projection bounds to
        // scaffold coordinate system
        int left_projection = -1;
        int right_projection = -1;
        int orient = -1;
        int iterator = 0;
        orient = get_orient(aln->contig_id, aln->strand, start_gap);
        if (aln->contig_id == cid_1) {
            // the contig to the left of the gap is aligned to the projector
            iterator = 1;
            if (aln->strand == '+') {
                int contig_origin = sgap->start - aln->s_length + 1;
                left_projection = contig_origin + near_projection;
                right_projection = contig_origin + far_projection;
            } else {
                int contig_origin = sgap->start;
                left_projection = contig_origin - near_projection;
                right_projection = contig_origin - far_projection;
            }
        } else if (aln->contig_id == cid_2) {
            // the contig to the right of the gap is aligned to the projector
            iterator = -1;
            if (aln->strand == '+') {
                int contig_origin = sgap->end + aln->s_length - 1;
                left_projection = contig_origin - far_projection;
                right_projection = contig_origin - near_projection;
            } else {
                int contig_origin = sgap->end;
                left_projection = contig_origin + far_projection;
                right_projection = contig_origin + near_projection;
            }
        } else {
            DIE("should never get here\n");
        }
        // The orientation of the projected read is opposite the projector
        orient = (orient == '+') ? '-' : '+';
        // Iterate over gaps, projecting mate into each gap that is within range
        // Note that we can be iterating up or down, i.e. iterator can be 1 or -1
        // Each thread is allocated a fixed block of gaps of size _ngaps_block, but only 
        // _gaps[MYTHREAD].ngaps of those are used. So to iterate over all gaps, we need to skip 
        // over the unused entries at the end of each thread block. We do this setting the
        // get_gap function to return NULL when a gap is out of range of allocated gaps
        for (; test_gapi < _max_ngaps && test_gapi >= 0; test_gapi += iterator) {
            shared [] gap_t *gap = get_gap(test_gapi, __LINE__, NULL, NULL_ON_FAIL);
            // only process if we have are in the range of gaps for this thread
            if (gap) {
                if (gap->scaff_id != scaff_id)
                    break;
                int placed = 0;
                int left = gap->start;
                int right = gap->end;
                if (gap->start >= gap->end) {
                    left = gap->end;
                    right = gap->start;
                }
                if (right_projection > left && left_projection < right) 
                    placed = 1;
                else if ((iterator == 1 && left > right_projection) || 
                         (iterator == -1 && right < left_projection)) 
                    break;
                if (placed) {
                    add_to_read_list(read_list, mate_name, test_gapi, orient);
                    nprojected++;
                }
            }
        }
        free(mate_name);
    }
    free(iter);
    destroy_htable(reads_to_alns, FREE_VALS);
    stop_timer(T_PROJECTED_ALNS);
    return nprojected;
}


static int get_pair_aln(pair_alns_t *pair_alns, char *buf, long *rejected, 
                        int reverse_complement, int five_prime_wiggle, int three_prime_wiggle)
{
    alignment_t aln;

    int res = get_alignment(buf, &aln, pair_alns->name);

    if (res != GOOD_ALIGNMENT) {
        rejected[res]++;
        return pair_alns->nalns;
    }
    // Assess alignment for completeness (do this before scaffold 
    // coordinate conversion!)
    // Important: Truncations are performed before reverse complementation
    // and apply to the end of the actual read
    int unaligned_start = aln.q_start - 1;
    int projected_start = (aln.strand == '+' ? 
                           aln.s_start - unaligned_start :
                           aln.s_stop + unaligned_start);
    aln.start_status = classify_alignment(unaligned_start, projected_start, aln.s_length,
                                          five_prime_wiggle, _cfg->truncate, 5);
    if (!aln.start_status) {
        rejected[REJECT_5_TRUNCATED]++;
        return pair_alns->nalns;
    }

    int unaligned_end = aln.q_length - aln.q_stop;
    int projected_end = (aln.strand == '+' ? 
                         aln.s_stop + unaligned_end :
                         aln.s_start - unaligned_end);
    aln.end_status = classify_alignment(unaligned_end, projected_end, aln.s_length, 
                                        three_prime_wiggle, _cfg->truncate, 3);
    if (!aln.end_status) {
        rejected[REJECT_3_TRUNCATED]++;
        return pair_alns->nalns;
    }
    // Re-orient alignment if requested
    if (reverse_complement) {
        aln.strand = (aln.strand == '+' ? '-' : '+');
        aln.q_start = aln.q_length - aln.q_stop + 1;
        aln.q_stop = aln.q_length - aln.q_start + 1;
    }
    pair_alns->alns[pair_alns->nalns] = aln;
    pair_alns->nalns++;
    return pair_alns->nalns;
}    

// returns the base name of a pair: @pair from any: '@pair/1' '@pair 1:Y:...' etc
static int get_base_name(char *name, char *buf)
{
    char *start_name = strchr(buf, '@');
    if (!start_name)
        return 0;
    assert(start_name != buf);
    assert(*(start_name - 1) == '\t');

    char *end_name = strchr(start_name, '\t');
    if (!end_name) 
        DIE("Invalid format for alignment: %s\n", buf);
    int len = (end_name - start_name) - 1;
    if (len > 3 && start_name[len - 2] == '/' && (start_name[len-1] == '1' || start_name[len-1] == '2')) {
        // @readpair/1 or @readpair/2  return @readpair
        len -= 2;
    }
    strncpy(name, start_name, len);
    name[len] = '\0';
    return 1;
}

static int read_pair_alns(FILE *f, char *fname, pair_alns_t *pair_alns, char *buf, char *curr, 
                          long *line, long *rejected, int reverse_complement, 
                          int five_prime_wiggle, int three_prime_wiggle, long *excessive_alns)
{
    if (feof(f))
        return 0;
    // this happens on the first invocation only
    if (!buf[0]) {
        if (!fgets(buf, MAX_BM_LINE_LEN - 1, f)) 
            return 0;
    }
    char name[MAX_READ_NAME_LEN] = "", curr_name[MAX_READ_NAME_LEN] = "";
    do {
        CHECK_LINE_LEN(buf, MAX_BM_LINE_LEN, fname);
        (*line)++;
        if (!get_base_name(name, buf)) 
            DIE("Could not find %s at line %ld in %s\n", name, *line, fname);
        if (curr[0]) {
            if (!get_base_name(curr_name, curr))
                DIE("Could not find %s in %s in %s\n", curr_name, curr, fname);
            if (strcmp(curr_name, name) == 0) {
                // accumulate alns
                if (pair_alns->nalns == MAX_PAIR_ALNS) 
                    (*excessive_alns)++;
                else 
                    get_pair_aln(pair_alns, buf, rejected, reverse_complement, 
                                 five_prime_wiggle, three_prime_wiggle);
            } else {
                // switching to new read
                strcpy(curr, buf);
                return 1;
            }
        } else {
            strcpy(curr, buf);
            get_pair_aln(pair_alns, buf, rejected, reverse_complement, 
                         five_prime_wiggle, three_prime_wiggle);
        }
    } while (fgets(buf, MAX_BM_LINE_LEN-1, f));
    return 1;
}


#define MIN(a, b) ((a) < (b) ? (a) : (b))

// Reads in alignments. Each line consists of:
// BLAST_TYPE  QUERY Q_START Q_STOP Q_LENGTH SUBJECT S_START S_STOP S_LENGTH STRAND SCORE E_VALUE IDENTITIES ALIGN_LENGTH
// The subject is the contig id.
// Expects the file to be split, one per thread
void process_meraligner_file(char *library_name, int reverse_complement, int insert_size, 
                             int insert_sigma, int five_prime_wiggle, int three_prime_wiggle)
{
    serial_printf("Processing %s/%s-%s_*_Read*... ", _cfg->base_dir, library_name, _cfg->mer_file);
#ifdef STEAL_PLACEMENTS
    serial_printf("Using work stealing to balance placements\n");
#endif
    START_TIMER(T_PROCESS_MERALIGNER);
    // guess 50 reads per gap - will be dynamically increased if needed
    read_list_t read_list = { .num = 0, .size = 50 * _ngaps_block };
    read_list.reads = calloc_chk(read_list.size, sizeof(read_orient_t));

    long *rejected = calloc_chk(NUM_REJECT_REASONS, sizeof(long));
    long num_aligned_gap_reads = 0;
    long num_projected_gap_reads = 0;
    long nalns = 0;
    // expect two files per library per thread, with endings alternating Read1 and Read2
    char fname[2][MAX_FILE_PATH];
    FILE *f[2] = {NULL, NULL};
    for (int i = 0; i < 2; i++) {
        sprintf(fname[i], "%s/%s-%s_%d_Read%d", 
                _cfg->base_dir, library_name, _cfg->mer_file, MYTHREAD, i+1);
        get_rank_path(fname[i], MYTHREAD);
        if (!i) 
            f[i] = fopen_chk(fname[i], "r");
        else 
            // unpaired reads can still be used too but will not have a Read1 file
            f[i] = fopen(fname[i], "r");
    }
    if (!f[1])
        tprintf("Could not open second file %s... using unpaired reads\n", fname[1]);

    tprintf_flush("Processing %s (+1)... ", fname[0]);

    long line0 = 0, line1 = 0;
    char buf0[MAX_BM_LINE_LEN] = "", buf1[MAX_BM_LINE_LEN] = "";
    char curr0[MAX_BM_LINE_LEN] = "", curr1[MAX_BM_LINE_LEN] = "";

    // count lines in file to determine maximum alignment blocks needed
    long max_pair_alns = 0;
    while (fgets(buf0, MAX_BM_LINE_LEN - 1, f[0])) 
        max_pair_alns++;
    rewind(f[0]);
    buf0[0] = 0;

    shared pair_alns_array_t *pair_alns_array = upc_all_alloc(THREADS, sizeof(pair_alns_array_t));
#ifdef STEAL_PLACEMENTS
    // there is one finished thread count per thread to avoid all threads hitting the same shared
    // variable when checking for termination
    shared int *finished_threads = upc_all_alloc(THREADS, sizeof(int));
    finished_threads[MYTHREAD] = 0;
#endif
    upc_barrier;
    pair_alns_array[MYTHREAD].pairs = upc_alloc(sizeof(pair_alns_t) * max_pair_alns);
    pair_alns_array[MYTHREAD].npairs = 0;
    pair_alns_array[MYTHREAD].next_pair = 0;

    int npairs;
    int max_nalns = 0;
    long excessive_alns = 0;
    for (npairs = 0; npairs < max_pair_alns; npairs++) {
        pair_alns_t *pair_alns = (pair_alns_t*)&pair_alns_array[MYTHREAD].pairs[npairs];
        pair_alns->nalns = 0;
        pair_alns->name[0] = 0;
        if (!read_pair_alns(f[0], fname[0], pair_alns, buf0, curr0, &line0, rejected, 
                            reverse_complement, five_prime_wiggle, three_prime_wiggle, 
                            &excessive_alns))
            break;
        if (f[1] && !read_pair_alns(f[1], fname[1], pair_alns, buf1, curr1, &line1, 
                                    rejected, reverse_complement, five_prime_wiggle, 
                                    three_prime_wiggle, &excessive_alns))
            DIE("EOF before finding matching pair for %s\n", pair_alns->name);
        if (max_nalns < pair_alns->nalns)
            max_nalns = pair_alns->nalns;
    }
    fclose(f[0]);
    if (f[1] != NULL) 
        fclose(f[1]);
    pair_alns_array[MYTHREAD].npairs = npairs;
    tprintf_flush("... done reading %ld aln pairs (max %ld )\n", npairs, max_pair_alns);
    tprintf_flush("Maximum alignments per pair: %d\n", max_nalns);
    double gb_to_alloc = (double)(sizeof(pair_alns_t) * max_pair_alns) / ONE_GB;
    tprintf_flush("Allocated %.3f GB for %ld pairs\n", gb_to_alloc, max_pair_alns);
    if (excessive_alns)
        tprintf_flush("Pairs with alignments > %d: %ld\n", MAX_PAIR_ALNS, excessive_alns);

    // we need a barrier here because we could steal from other threads
    upc_barrier;

    // The stealing strategy is to randomly pick a victim initially to steal from, and then
    // iterate from that victim onwards for subsequent steals. Also, for each steal, 20 victims
    // are checked, and the one with the most available to steal is chosen.
    // To make sure we terminate, each thread increments a counter to indicate that it has 
    // finished all its local entries, so it knows when to break out of the loop

    unsigned rseed = MYTHREAD;
    srand(rseed);
    int stealing = 0;
    pair_alns_t *pair_alns;
    pair_alns_t pair_alns_block[STEAL_BLOCK];
    long next_pair;
    long num_pairs_in_block;
    long num_steals = 0;
    long tot_num_pairs = 0;
    while (1) {
        if (!stealing) {
            next_pair = bupc_atomicL_fetchadd_relaxed(&(pair_alns_array[MYTHREAD].next_pair), STEAL_BLOCK);
            if (next_pair >= pair_alns_array[MYTHREAD].npairs) {
#ifdef STEAL_PLACEMENTS
                // update all finished thread counts on all other threads
                for (int i = 0; i < THREADS; i++) 
                    bupc_atomicI_fetchadd_relaxed(&finished_threads[i], 1);
                stealing = 1;
#else
                break;
#endif
                
            } else {
                num_pairs_in_block = MIN(pair_alns_array[MYTHREAD].npairs - next_pair, STEAL_BLOCK);
            }
        }
#ifdef STEAL_PLACEMENTS
        if (stealing) {
            if (bupc_atomicI_read_relaxed(&finished_threads[MYTHREAD]) == THREADS)
                break;
            num_pairs_in_block = 0;
            int victim = -1, max_remainder = 0;
            // find slowest target thread in random subset
            for (int i = 0; i < 20; i++) {
                int new_victim = rand_r(&rseed) % THREADS;
                if (new_victim != MYTHREAD) {
                    int remainder = pair_alns_array[new_victim].npairs - 
                        bupc_atomicL_read_relaxed(&(pair_alns_array[new_victim].next_pair));
                    if (remainder > max_remainder) {
                        max_remainder = remainder;
                        victim = new_victim;
                    }
                }
            }
            if (victim != -1) {
                next_pair = bupc_atomicL_fetchadd_relaxed(&(pair_alns_array[victim].next_pair), STEAL_BLOCK);
                if (next_pair < pair_alns_array[victim].npairs - STEAL_BLOCK) {
                    num_pairs_in_block = MIN(pair_alns_array[victim].npairs - next_pair, STEAL_BLOCK);
                    // copy the remote data to local
                    upc_memget(&pair_alns_block, &(pair_alns_array[victim].pairs[next_pair]), 
                               sizeof(pair_alns_t) * num_pairs_in_block);
                    num_steals += num_pairs_in_block;
                }
            }                
        }
#endif
        tot_num_pairs += num_pairs_in_block;
        for (int pi = 0; pi < num_pairs_in_block; pi++) {
            if (!stealing)
                pair_alns = (pair_alns_t*)&(pair_alns_array[MYTHREAD].pairs[pi + next_pair]);
            else
                pair_alns = &(pair_alns_block[pi]);
            if (pair_alns->nalns) {
                num_aligned_gap_reads += direct_alignments(pair_alns, &read_list);
                if (_cfg->pair_projection && f[1] > 0)
                    num_projected_gap_reads += projected_alignments(pair_alns, &read_list, 
                                                                    insert_size, insert_sigma);
                nalns += pair_alns->nalns;
            }
        }
    }
    upc_barrier;
    upc_free(pair_alns_array[MYTHREAD].pairs);
    upc_barrier;
    upc_all_free(pair_alns_array);
#ifdef STEAL_PLACEMENTS
    upc_all_free(finished_threads);
#endif

    long tot_unused = 0;
    for (int i = 0; i < NUM_REJECT_REASONS; i++) 
        tot_unused += rejected[i];
    tprintf_flush("Found %ld alignments found, used %.3f fraction\n", 
				  nalns, (double)nalns / (nalns + tot_unused));
    tprintf_flush("Number of steals: %ld (%.1f %%)\n", 
                  num_steals, (100.0 * num_steals / tot_num_pairs));
    tprintf_flush("Unused alignments:\n");
    for (int i = 0; i < NUM_REJECT_REASONS; i++) 
        if (rejected[i]) 
            tprintf_flush("\t%s\t%d\n", reject_reasons[i], rejected[i]);
    tprintf_flush("Total reads placed in gaps = %d (aligned) + %d (projected)\n", 
                  num_aligned_gap_reads, num_projected_gap_reads);
    free(rejected);
    tprintf_flush("Direct placement took %.3f s, projected placement took %.3f s\n", 
                  get_elapsed_time(T_DIRECT_ALNS), get_elapsed_time(T_PROJECTED_ALNS));

    stop_timer(T_PROCESS_MERALIGNER);
    tprintf_flush("Placing reads took %.3f s\n", get_elapsed_time(T_PROCESS_MERALIGNER));
    START_TIMER(T_PROCESS_MERALIGNER);

    // now put reads in shared hash table
    long tot_reads = reduce_long(read_list.num, UPC_ADD);
    tprintf_flush("Found %ld reads, total %ld\n", read_list.num, tot_reads);

    dhtable_init(tot_reads, READ_HASH_TABLE_LOAD_FACTOR);

    START_TIMER(T_PUT_READ_ORIENT);
    tprintf_flush("Number of reads in list %d and number of reads allocated %d (reallocs %d)\n",
                  read_list.num, read_list.size, _num_reallocs);
    long num_put_failures = 0;
    long num_puts = 0;
    char name_buf[MAX_READ_NAME_LEN + 2];
    for (long i = 0; i < read_list.num; i++) {
        num_puts++;
        read_orient_t *read_orient = &read_list.reads[i];
        sprintf(name_buf, "%c%s", read_orient->orient, read_orient->name);
        if (!dhtable_put(name_buf, read_orient->gap_id))
            num_put_failures++;
    }
	upc_synci();
    free(read_list.reads);

    long tot_num_put_failures = reduce_long(num_put_failures, UPC_ADD);
    if (tot_num_put_failures) {
        long tot_num_puts = reduce_long(num_puts, UPC_ADD);
        double perc_fails = (100.0 * tot_num_put_failures) / tot_num_puts;
        if (!MYTHREAD && perc_fails >= 5) {
            serial_printf("\n");
            WARN("\nDistributed hash table had %ld put failures out of %ld puts (%.2f %%).\n"
                 "[Indicates potential loss of efficiency and possibly fewer gaps closed]\n", 
                 tot_num_put_failures, tot_num_puts, perc_fails);
        }
        //check_read_data_hash();
    }
    stop_timer(T_PUT_READ_ORIENT);

    stop_timer(T_PROCESS_MERALIGNER);
    print_memory_used();
    serial_printf("done in %.3f s\n", get_elapsed_time(T_PROCESS_MERALIGNER));
    UPC_TIMED_BARRIER;
    //double t_all = get_timer_all_max(T_PROCESS_MERALIGNER);
    //serial_printf("All done in %.3f s\n", t_all);
#ifdef CHECK_READ_DATA_HASH
    check_read_data_hash();
#endif
}

char *get_primer(shared contig_t *contig, int contig_ext, int reverse_flag)
{
    char *primer = malloc_chk(_cfg->mer_size + 1);
    if (contig_ext == 5) {
        upc_memget(primer, contig->seq, _cfg->mer_size);
        primer[_cfg->mer_size] = 0;
    } else {
        char *seq = malloc_chk(contig->seq_len + 1);
        UPC_MEMGET_STR(seq, contig->seq, contig->seq_len + 1);
        int offset = contig->seq_len - _cfg->mer_size;
        strcpy(primer, seq + offset);
        primer[_cfg->mer_size] = 0;
        free(seq);
    }
    if (contig_ext == reverse_flag)
        switch_code(reverse(primer));
    return primer;
}

void print_gaps(void) 
{
    serial_printf("Printing gaps... ");
    START_TIMER(T_PRINT_GAPS);
    char fname[MAX_FILE_PATH];
    sprintf(fname, "%s/gaps.%d", _cfg->base_dir, MYTHREAD);
    get_rank_path(fname, MYTHREAD);
    FILE *f = fopen_chk(fname, "w+");
    long max_nreads = 0;
    long tot_nreads = 0;
    for (int i = 0; i < _gaps[MYTHREAD].ngaps; i++) {
        gap_t *gap = (gap_t*)&(_gaps[MYTHREAD].gaps[i]);
        tot_nreads += gap->nreads;
        if (gap->nreads > max_nreads)
            max_nreads = gap->nreads;
        shared contig_t *contig1;
        get_contig_chk(contig1, gap->contig1_id);
        shared contig_t *contig2;
        get_contig_chk(contig2, gap->contig2_id);
        char *primer1 = get_primer(contig1, gap->contig_ext1, 5);
		if (!primer1)
			DIE("Could not find primer1 for gap %d\n", i);
        char *primer2 = get_primer(contig2, gap->contig_ext2, 3);
		if (!primer2)
			DIE("Could not find primer2 for gap %d\n", i);

        char buf[1000];
        sprintf(buf, "Scaffold%d\tContig%d.%d\t%s\tContig%d.%d\t%s\t%d\t%.0f\t%d\t",
                gap->scaff_id, gap->contig1_id, gap->contig_ext1,
                primer1, gap->contig2_id, gap->contig_ext2,
                primer2, gap->size, gap->uncertainty, gap->nreads);
        for (int j = 0; j < gap->nreads; j++)
            fprintf(f, "%s %s:%s\n", buf, (char*)GET_READ(gap->reads_nts, i), 
                    (char*)GET_READ(gap->reads_quals, j));
        /*
		  fprintf(f, "Scaffold%d\tContig%d.%d\t%s\tContig%d.%d\t%s\t%d\t%.0f", 
		  gap->scaff_id, gap->contig1_id, gap->contig_ext1, 
		  primer1, gap->contig2_id, gap->contig_ext2, 
		  primer2, gap->size, gap->uncertainty);
		  for (int j = 0; j < gap->nreads; j++) 
		  fprintf(f, "\t%s:%s", gap->reads[j].nts, gap->reads[j].quals);
		  fprintf(f, "\n");
        */
		free(primer1);
		free(primer2);
    }
	fclose(f);
    stop_timer(T_PRINT_GAPS);
	serial_printf("done in %.3f s\n", get_elapsed_time(T_PRINT_GAPS));
}

static int get_reads_for_orient(char orient, char *read_name, int *gap_ids, char *curr_read_nts, 
                                char *curr_read_quals)
{
    char buf[MAX_READ_NAME_LEN + 2];
    sprintf(buf, "%c%s", orient, read_name);
    // the dhtable can have multiple enties per read, so we get them all in an array of gap ids
    int ngap_ids = dhtable_get(buf, gap_ids, MAX_READ_GAP_IDS);
    if (!ngap_ids)
        return 0;
    if (orient == '-') {
        switch_code(reverse(curr_read_nts));
        reverse(curr_read_quals);
    }
    int curr_read_len = strlen(curr_read_nts);
    START_TIMER(T_GAP_COPY);
    for (int i = 0; i < ngap_ids; i++) {
        int gi = gap_ids[i];
        shared [] gap_t *gap = get_gap(gi, __LINE__, NULL, ABORT_ON_FAIL);
        int max_reads = gap->max_reads;
        if (!max_reads) 
            continue;
        int ri = bupc_atomicI_fetchadd_relaxed(&(gap->nreads), 1);
        if (ri >= max_reads) {
            char buf[256];
            WARN("max reads exceeded for gap %d, %s, max is %d\n"
                 "Duplicate reads in FASTQ? Latest read:\n%s\n", 
                 gi, get_gapi_info(gi, buf), max_reads, read_name);
            continue;
        }
        if (!gap->reads_lens || !gap->reads_nts || !gap->reads_quals)
            DIE("gap->reads is NULL for gap %d\n", gi);
        gap->reads_lens[ri] = curr_read_len;
        upc_memput_nbi(GET_READ(gap->reads_nts, ri), curr_read_nts, curr_read_len + 1);
        upc_memput_nbi(GET_READ(gap->reads_quals, ri), curr_read_quals, curr_read_len + 1);
    }
    stop_timer(T_GAP_COPY);
    return 1;
}

void process_fastq_file(char *library_name)
{
    START_TIMER(T_PROCESS_FASTQ);
    // allocate local shared for the reads for each gap
    int num_no_read_gaps = 0;
    for (int i = 0; i < _gaps[MYTHREAD].ngaps; i++) {
        gap_t *gap = (gap_t*)&(_gaps[MYTHREAD].gaps[i]);

        int max_reads = gap->max_reads;
        if (max_reads > MAX_READS_PER_GAP) {
            tprintf_flush("WARNING: excluding gap %d: too many reads %d > %d\n", 
                          gap->id, max_reads, MAX_READS_PER_GAP);
            gap->max_reads = 0;
            continue;
        }
        if (max_reads == 0) {
            num_no_read_gaps++;
            continue;
        }
        if (max_reads < 0)
            DIE("negative max_reads %d\n", max_reads);

#ifdef DEBUG
        tprintf("For library %s, gap %d, allocating %ld reads of %d max_read_len\n", library_name, gap->id, max_reads, _cfg->max_read_len);
#endif

        // FIXME: this is inefficient: since there is no realloc in upc, we allocate a new 
        // array, copy across the old one, and then delete the old one. Ugh.
        shared [] int *new_reads_lens = upc_alloc(max_reads * sizeof(int));
        shared [] char *new_reads_nts = upc_alloc(max_reads * (_cfg->max_read_len + 1));
        shared [] char *new_reads_quals = upc_alloc(max_reads * (_cfg->max_read_len + 1));
        if (gap->nreads) {
#ifdef DEBUG
            tprintf("  copying across %ld old reads\n", gap->nreads);
#endif
            // copy across old reads
            assert(gap->reads_lens && IS_VALID_UPC_PTR(gap->reads_lens));
            assert(gap->reads_nts && IS_VALID_UPC_PTR(gap->reads_nts));
            assert(gap->reads_quals && IS_VALID_UPC_PTR(gap->reads_quals));
            memcpy((char*)new_reads_lens, (char*)gap->reads_lens, gap->nreads * sizeof(int));
            memcpy((char*)new_reads_nts, (char*)gap->reads_nts, gap->nreads * (_cfg->max_read_len + 1));
            memcpy((char*)new_reads_quals, (char*)gap->reads_quals, gap->nreads * (_cfg->max_read_len + 1));
            upc_free(gap->reads_lens); gap->reads_lens = NULL;
            upc_free(gap->reads_nts); gap->reads_nts = NULL;
            upc_free(gap->reads_quals); gap->reads_quals = NULL;
        }
        gap->reads_lens = new_reads_lens;
        gap->reads_nts = new_reads_nts;
        gap->reads_quals = new_reads_quals;
    }
    UPC_TIMED_BARRIER;
    if (num_no_read_gaps)
        tprintf("Number of gaps without reads (found in srf but not meraligner): %d (%.2f fraction)\n", 
                num_no_read_gaps, (double)num_no_read_gaps / _gaps[MYTHREAD].ngaps);

    long num_reads = 0, num_reads_used = 0;
    int gap_ids[MAX_READ_GAP_IDS];
    char *read_name = NULL, *curr_read_nts = NULL, *curr_read_quals = NULL;
    Buffer idBuf = initBuffer(MAX_READ_NAME_LEN), seqBuf = initBuffer(MAX_READ_LEN), 
        qualBuf = initBuffer(MAX_READ_LEN);
    int curr_read_len;

    int from_shm = 0;
    if (strstr(_cfg->base_dir, "/dev/shm"))
        from_shm = 1;

    // Get FileOfFileNames by library 
    char fofn[MAX_FILE_PATH], fname[MAX_FILE_PATH];
    snprintf(fofn, MAX_FILE_PATH, "%s.fofn", library_name);
    serial_printf("Processing fastq files from %s...\n", fofn);
    tprintf_flush("Processing fastq files from %s...\n", fofn);
    FILE *fofnFD = fopen_chk(fofn, "r");
    while (fgets(fname, MAX_FILE_PATH, fofnFD)) { 
        fname[strlen(fname)-1] = '\0';

        fq_reader_t fqr = create_fq_reader();
        open_fq(fqr, fname, from_shm);
        
        while (get_next_fq_record(fqr, idBuf, seqBuf, qualBuf)) {
            if (getLengthBuffer(seqBuf) > _cfg->max_read_len)
                DIE("Line length in FASTQ file is higher than maximum from -l parameter: %ld > %d\n", 
                    getLengthBuffer(seqBuf), _cfg->max_read_len);
            read_name = getStartBuffer(idBuf);
            curr_read_nts = getStartBuffer(seqBuf); 
            curr_read_quals = getStartBuffer(qualBuf);
            num_reads++;

#ifdef DBG_FASTQ
            if (read_name[strlen(read_name) - 1] == '2') {
                dbg("%s\t", read_name);
                dbg("%s\t", curr_read_nts);
                dbg("+%s\t", read_name + 1);
                //dbg("+\n");
                dbg("%s\n", curr_read_quals);
            }
#endif
            num_reads++;
            // try both orientations
            int use_read = 0;
            if (get_reads_for_orient('+', read_name, gap_ids, curr_read_nts, curr_read_quals))
                use_read = 1;
            if (get_reads_for_orient('-', read_name, gap_ids, curr_read_nts, curr_read_quals))
                use_read = 1;
            num_reads_used += use_read;
        } // reading fastq 
        destroy_fq_reader(fqr);
    } // reading fofn 

    // sync for all the nonblocking gap read puts
    upc_synci();

    freeBuffer(seqBuf);
    freeBuffer(qualBuf);

    UPC_TIMED_BARRIER;

    dhtable_free();
    
    long num_missing_reads = 0;
    for (int i = 0; i < _gaps[MYTHREAD].ngaps; i++) {
        gap_t *gap = (gap_t*)&(_gaps[MYTHREAD].gaps[i]);
        if (gap->max_reads < gap->nreads)
            gap->nreads = gap->max_reads;
        num_missing_reads += (gap->max_reads - gap->nreads);
#ifdef DBG_DREAD_AFTER
        char buf[200];
        for (int j = 0; j < gap->nreads; j++) 
            dbg("DREADAFTER %s %s\n", get_gapi_info(gap->id, buf), (char*)gap->reads_nts[j]);
#endif            
    }
    if (num_missing_reads) {
        double perc_missing = (100.0 * num_missing_reads) / num_reads;
        if (perc_missing > 5)
            tprintf_flush("WARNING: found %ld (%.2f %%) reads in meraligner files that are either "
                          "not in FASTQ files or not in hash table\n", 
                          num_missing_reads, perc_missing);
    }

    tprintf_flush("In FASTQ: %d reads found, used %.3f fraction, in %.3f s\n", 
                  num_reads, (double)num_reads_used / num_reads, get_elapsed_time(T_PROCESS_FASTQ));
    stop_timer(T_PROCESS_FASTQ);
    print_memory_used();
    UPC_TIMED_BARRIER;
    //double t_all = get_timer_all_max(T_PROCESS_FASTQ);
    serial_printf("Done with FASTQ files in %.3f s\n", get_elapsed_time(T_PROCESS_FASTQ));
}

static void correct_contig_gapi(int contig_id, int old_gi, int new_gi)
{
    shared contig_t *contig;
    get_contig_chk(contig, contig_id);
    if (!contig) 
        DIE("Missing contig %d for gap %d\n", contig_id, old_gi);
    if (contig->gapi_5 == old_gi)
        contig->gapi_5 = new_gi;
    else if (contig->gapi_3 == old_gi)
        contig->gapi_3 = new_gi;
    else 
        DIE("Could not find gap id %d (to replace with %d) for contig %d (%d),"
            " gapi_5 %d gapi_3 %d\n", 
            old_gi, new_gi, contig->id, contig_id, contig->gapi_5, contig->gapi_3);
}


void load_balance_gaps(void)
{
    START_TIMER(T_BALANCE_GAPS);
    serial_printf("Load balancing gaps...\n");
    
    shared int *steal_flags = upc_all_alloc(THREADS, sizeof(int));
    shared int *req_ngaps = upc_all_alloc(THREADS, sizeof(int));
    req_ngaps[MYTHREAD] = _tot_ngaps / THREADS;
    // ensure that the required gaps per thread does not differ by more than one
    if (MYTHREAD < _tot_ngaps - req_ngaps[MYTHREAD] * THREADS)
        req_ngaps[MYTHREAD]++;
    if (_gaps[MYTHREAD].ngaps <= req_ngaps[MYTHREAD])
        steal_flags[MYTHREAD] = 1;
    else 
        steal_flags[MYTHREAD] = 0;
    //int steal_block = 10;
    UPC_TIMED_BARRIER;
    int under_ngaps = req_ngaps[MYTHREAD] - _gaps[MYTHREAD].ngaps;
    int ngaps_before = _gaps[MYTHREAD].ngaps;
    int i = MYTHREAD + 1;
    while (under_ngaps > 0) {
        if (i == MYTHREAD)
            i++;
        if (i == THREADS)
            i = 0;
        if (!bupc_atomicI_cswap_relaxed(&steal_flags[i], 0, 1)) {
            // now steal
            int steal_num = _gaps[i].ngaps - req_ngaps[i];
            if (steal_num > under_ngaps)
                steal_num = under_ngaps;
            // limit the maximum number stolen at any time in order to better distribute the gaps
            //if (steal_num > steal_block)
            //    steal_num = steal_block;
            _gaps[i].ngaps -= steal_num;
            _gaps[i].gap_end = _gaps[i].gap_start + _gaps[i].ngaps;
            int steal_from_index = _gaps[i].ngaps;
            if (_gaps[i].ngaps > req_ngaps[i]) {
                // set flag to 0 to indicate that there's still some excess gaps at victim
                if (!bupc_atomicI_cswap_relaxed(&steal_flags[i], 1, 0))
                    DIE("Flag was set to 0 when it should have been set to 1 for stealing from %d\n", i);
            }
            // we've released the "lock" (cswap flag), so now copy across the data
            // Make sure to set the gap ids and contig ids correctly
            for (int gi = 0; gi < steal_num; gi++) {
                int srci = gi + steal_from_index;
                int dsti = gi + _gaps[MYTHREAD].ngaps;
                upc_memcpy(&_gaps[MYTHREAD].gaps[dsti], &_gaps[i].gaps[srci], sizeof(gap_t));
                gap_t *gap = (gap_t*)&(_gaps[MYTHREAD].gaps[dsti]);
                int new_gi = dsti + _gaps[MYTHREAD].gap_start;
                correct_contig_gapi(gap->contig1_id, gap->id, new_gi);
                correct_contig_gapi(gap->contig2_id, gap->id, new_gi);
                gap->id = new_gi;
            } 

            //int prev_ngaps = _gaps[MYTHREAD].ngaps;
            _gaps[MYTHREAD].ngaps += steal_num;
            _gaps[MYTHREAD].gap_end = _gaps[MYTHREAD].gap_start + _gaps[MYTHREAD].ngaps;
            under_ngaps = req_ngaps[MYTHREAD] - _gaps[MYTHREAD].ngaps;
            //tprintf_flush("Had %d gaps and Stole %d gaps from thread %d (with %d gaps)\n",
            //              prev_ngaps, steal_num, i, _gaps[i].ngaps);
        }
        i++;
    }
    UPC_TIMED_BARRIER;
    upc_all_free(steal_flags);
    upc_all_free(req_ngaps);
    tprintf_flush("Balanced gaps: %d -> %d\n", ngaps_before, _gaps[MYTHREAD].ngaps);
    stop_timer(T_BALANCE_GAPS);
    serial_printf("Done with gap balancing in %.3f s\n", get_elapsed_time(T_BALANCE_GAPS));
}
