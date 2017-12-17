#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <upc.h>
#include <upc_nb.h>
#include <upc_tick.h>

#include "utils.h"
#include "tracing.h"
#include "../common/optlist.h"

#define TICKS_TO_S(t) ((double)upc_ticks_to_ns(t) / 1000000000.0)
#define ELAPSED_TIME(t) TICKS_TO_S(upc_ticks_now() - (t))

typedef struct assembly_elem {
    int id;
    int len;
    char *seq;
    struct assembly_elem *next;
} assembly_t;

typedef struct {
    int id;
    int len;
    shared [] char *seq;
} shared_assembly_t;

typedef struct {
    shared [] shared_assembly_t *a;
    long num;
} assembly_set_t;

#define FASTA_LINE_WIDTH 100

static char _fa_fname[MAX_FILE_PATH] = "assembly.";
static char _fa_fname_suffix[MAX_FILE_PATH] = ".fa";
static char _output_fname[MAX_FILE_PATH] = "final_assembly.fa";

/*
static void print_fasta(FILE *f, int scaff_id, char *seq, int seq_len)
{
    fprintf(f, ">Scaffold%d\n", scaff_id);
    int nlines = INT_CEIL(seq_len, FASTA_LINE_WIDTH);
    char buf[FASTA_LINE_WIDTH + 1];
    for (int i = 0; i < nlines; i++) {
        strncpy(buf, seq + (i * FASTA_LINE_WIDTH), FASTA_LINE_WIDTH);
        buf[FASTA_LINE_WIDTH] = 0;
        fprintf(f, "%s\n", buf);
    }
}
*/

#define CMP_SEQS(len1, len2, s1, s2) ((len1) == (len2) ? strcmp((s2), (s1)) : (len2) - (len1))
#define CMP_ASSEMBLIES(a1, a2) CMP_SEQS((a1)->len, (a2)->len, (a1)->seq, (a2)->seq)

static int cmp_seqs(const void *p1, const void *p2)
{
    shared_assembly_t *a1 = (shared_assembly_t*)p1;
    shared_assembly_t *a2 = (shared_assembly_t*)p2;
    return CMP_SEQS(a1->len, a2->len, (char*)a1->seq, (char*)a2->seq);
}

static void revcomp(char *seq, int len)
{
    char *rev_seq = strdup(seq);
    switch_code(reverse(rev_seq));
    if (strcasecmp(rev_seq, seq) < 0) 
        strcpy(seq, rev_seq);
    free(rev_seq);
}

static int get_pos(assembly_t **assemblies, int start, int end, assembly_t *assembly)
{
    // binary search to find the position in the sequence
    int low = start;
    int high = end - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        assembly_t *midval = assemblies[mid];
        int cmp = CMP_ASSEMBLIES(midval, assembly);
        if (cmp < 0)
            low = mid + 1;
        else if (cmp > 0)
            high = mid - 1;
        else
            return mid; 
    }
    return low;
}

static assembly_t *sorted_merge(assembly_t *a, assembly_t *b)
{
    if (!a)
        return b;
    if (!b)
        return a;
    assembly_t dummy = { .next = NULL };
    assembly_t *tail = &dummy;
    while (a && b) {
        if (CMP_ASSEMBLIES(a, b) > 0) {
            tail->next = a;
            a = a->next;
        } else {
            tail->next = b;
            b = b->next;
        }
        tail = tail->next;
        tail->next = NULL;
    }
    if (!a)
        tail->next = b;
    if (!b)
        tail->next = a;
    return dummy.next;
}

static void merge_assemblies(shared assembly_set_t *assemblies, int cores_per_node) 
{
    // this is done only on thread 0
    long tot_anum = 0;
    upc_tick_t t_merge = upc_ticks_now();
    for (int ti = 0; ti < THREADS; ti += cores_per_node)
        tot_anum += assemblies[ti].num;
    assembly_t *all_assemblies = NULL;
    upc_tick_t t_fetch = 0;
    int num_seqs = 0;
    // merge all the data to local
    for (int ti = 0; ti < THREADS; ti += cores_per_node) {
        long remote_num = assemblies[ti].num;
        assembly_t *head = NULL;
        upc_tick_t t = upc_ticks_now();
        for (long i = 0; i < remote_num; i++) {
            assembly_t *assembly = malloc_chk(sizeof(assembly_t));
            assembly->len = assemblies[ti].a[i].len;
            assembly->id = assemblies[ti].a[i].id;
            assembly->seq = malloc_chk(assembly->len + 1);
            upc_memget_nbi(assembly->seq, assemblies[ti].a[i].seq, assembly->len + 1);
            assembly->next = head;
            head = assembly;
        }
        upc_synci();
        num_seqs += remote_num;
        t_fetch += (upc_ticks_now() - t);
        all_assemblies = sorted_merge(all_assemblies, head);
    }
    serial_printf("Merging took %.3f s\n", ELAPSED_TIME(t_merge) - TICKS_TO_S(t_fetch));
    serial_printf("Fetching took %.3f s\n", TICKS_TO_S(t_fetch));
    serial_printf("Found %d sequences\n", num_seqs);

    upc_tick_t t_write = upc_ticks_now();

    serial_printf("Writing canonical sequences to %s\n", _output_fname);
    FILE *f = fopen_chk(_output_fname, "w");
    int i = 0;
    long tot_size = 0;
    long *lens = malloc_chk(num_seqs * sizeof(long));
    for (assembly_t *cur = all_assemblies; cur; cur = cur->next) {
        fprintf(f, ">%d-%d\n%s\n", i, cur->len, cur->seq);
        tot_size += cur->len;
        lens[i] = cur->len;
        i++;
    }
    fclose(f);
    serial_printf("Writing took %.3f s\n", ELAPSED_TIME(t_write));
    // now compute N50 and L50. It's almost free since the scaffolds are already sorted
    long running_tot = 0;
    const int NUM_INDEXES = 5;
    double indexes[] = {0.1, 0.25, 0.5, 0.75, 0.9};
    int n_index = 0;
    serial_printf("Assembly stats:\n\tSize %ld\n", tot_size);
    for (i = num_seqs - 1; i >= 0 && n_index < NUM_INDEXES; i--) {
        running_tot += lens[i];
        if (running_tot >= indexes[n_index] * tot_size) {
            int nx = (int)(indexes[n_index] * 100);
            serial_printf("\tN%d %ld L%d %d\n", nx, lens[i], nx, i + 1);
            n_index++;
        }
    }
}

static void addCanonicalFasta(shared_assembly_t *assemblies, long *anum, char *seq) {
    assemblies[*anum].len = strlen(seq);
    if (assemblies[*anum].len <= 0)
        DIE("len is %d\n", assemblies[*anum].len);
    assemblies[*anum].seq = upc_alloc(assemblies[*anum].len + 1);
    revcomp(seq, assemblies[*anum].len);
    strcpy((char*)assemblies[*anum].seq, seq);
    // always make all upper case
    //for (int j = 0; j < assemblies[*anum].len + 1; j++)
    //    assemblies[*anum].seq[j] = toupper(seq[j]);
    (*anum)++;
}

static long count_assemblies(char *output_dir, int cores_per_node)
{
    char fname[MAX_FILE_PATH];
    char buf[FASTA_LINE_WIDTH + 2];
    long anum = 0;
    for (int i = MYTHREAD; i < MYTHREAD + cores_per_node; i++) {
        sprintf(fname, "%s/%s%d%s", output_dir, _fa_fname, i, _fa_fname_suffix);
        FILE *f = fopen_rank_path(fname, "r", i);
        char *res = NULL;
        do {
            res = fgets(buf, FASTA_LINE_WIDTH + 2, f);
            if (res == NULL) break;
            if (buf[0] =='>') 
                anum++;
        } while (res != NULL);
        fclose(f);
    }
    return anum;
}


static long read_assemblies(char *output_dir, shared_assembly_t *assemblies, 
                            long max_assemblies, int cores_per_node)
{
    for (long i = 0; i < max_assemblies; i++)
        assemblies[i].id = -1;
    char fname[MAX_FILE_PATH];
    char buf[FASTA_LINE_WIDTH + 2];
    int max_scaff_seq_len = 1000000;
    char *seq = malloc_chk(max_scaff_seq_len + 1);
    int seq_len = 0;
    int word_len = 0;
    long anum = 0;
    for (int i = MYTHREAD; i < MYTHREAD + cores_per_node; i++) {
        sprintf(fname, "%s/%s%d%s", output_dir, _fa_fname, i, _fa_fname_suffix);
        FILE *f = fopen_rank_path(fname, "r", i);
        char *res = NULL;
        do {
            res = fgets(buf, FASTA_LINE_WIDTH + 2, f);
            if (res == NULL) break;
            if (buf[0] == '>') {
                //printf("adding scaffold %ld\n", anum);
                if (!word_len)
                    word_len = (buf[1] == 'C' ? strlen(">Contig_") : strlen(">Scaffold"));
                //tprintf("read_assemblies: anum: %ld line: %s", anum, buf);
                if (assemblies[anum].id != -1) {
                    addCanonicalFasta(assemblies, &anum, seq);
                    if (anum > max_assemblies)
                        DIE("anum %ld > %ld\n", anum, max_assemblies);
                }
                if (res) {
                    seq[0] = 0;
                    seq_len = 0;
                    assemblies[anum].id = atoi(buf + word_len);
                }
            } else {
                if (assemblies[anum].id == -1)
                    DIE("corrupted assembly file, %s, got sequence before id\n", fname);
                int buf_len = strlen(buf);
                // get rid of newline
                if (buf[buf_len - 1] == '\n') {
                    buf_len--;
                    buf[buf_len] = 0;
                }
                int len = seq_len + buf_len;
                if (len >= max_scaff_seq_len) {
                    max_scaff_seq_len = 2 * len;
                    seq = realloc_chk(seq, max_scaff_seq_len + 1);
                }
                strcpy(seq + seq_len, buf);
                seq_len += buf_len;
            }
        } while (res != NULL);
        if (anum == max_assemblies)
            break;
        if (assemblies[anum].id != -1) 
            addCanonicalFasta(assemblies, &anum, seq);
        fclose(f);
    }
    free(seq);
    return anum;
}

int main(int argc, char **argv)
{
    upc_tick_t t_start = upc_ticks_now();
    char *output_dir = NULL;
    int max_seqs = 0;
    int cores_per_node = 0;
    option_t *thisOpt;
    option_t *optList = GetOptList(argc, argv, "o:O:s:n:f:F:");
    while (optList) {
        thisOpt = optList;
        optList = optList->next;
        switch (thisOpt->option) {
        case 'o':
            output_dir = thisOpt->argument;
            break;
        case 'O':
            strcpy(_output_fname, thisOpt->argument);
            break;
        case 's':
            max_seqs = atoi(thisOpt->argument);
            break;
        case 'n':
            cores_per_node = atoi(thisOpt->argument);
            break;
        case 'f':
            strcpy(_fa_fname, thisOpt->argument);
            break;
        case 'F':
            strcpy(_fa_fname_suffix, thisOpt->argument);
            break;
        default:
            serial_printf("Invalid option: %c\n", thisOpt->option);
            upc_barrier;
            upc_global_exit(1);
        }
    }
    if (!output_dir || !max_seqs || ! cores_per_node) {
        serial_printf("Usage: %s -o outputdir -s max_seqs -n cores_per_node -f fafiles\n",
                      argv[0]);
        upc_barrier;
        upc_global_exit(1);
    }
    set_output_dir(output_dir);
    serial_printf("Reading from %s/%s{0..%d}%s\n", output_dir, _fa_fname, THREADS-1, _fa_fname_suffix);
    set_log_prefix(_output_fname);
    shared assembly_set_t *assemblies = upc_all_alloc(THREADS, sizeof(assembly_set_t));
    upc_barrier;
    if (!(MYTHREAD % cores_per_node)) {
        long anum = count_assemblies(output_dir, cores_per_node);
        assemblies[MYTHREAD].a = upc_alloc(anum * sizeof(shared_assembly_t));
        upc_tick_t t_read = upc_ticks_now();
        read_assemblies(output_dir, (shared_assembly_t*)assemblies[MYTHREAD].a, 
                        anum, cores_per_node);
        assemblies[MYTHREAD].num = anum;
        serial_printf("Reading took %.3f s\n", ELAPSED_TIME(t_read));
        qsort((shared_assembly_t*)assemblies[MYTHREAD].a, anum, sizeof(shared_assembly_t), 
              cmp_seqs);
    }
    upc_barrier;
    serial_printf("Reading and qsort took %.3f s\n", ELAPSED_TIME(t_start));
    if (!MYTHREAD) 
        merge_assemblies(assemblies, cores_per_node);
    double elapsed_secs = ELAPSED_TIME(t_start);
    serial_printf("Canonical assembly took %.3f s\n", elapsed_secs);
    char buf[MAX_FILE_PATH];
    serial_printf("Overall time for %s is %.2f s\n", get_basename(buf, argv[0]), elapsed_secs);
    upc_barrier;
}
