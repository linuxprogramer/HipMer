/*
  Distributed hash table

 This form of hash table uses atomics and not locks. But it will only work if the gets and puts
 are not interleaved. This is the case here, as all the puts are done when reading the meraligner
 file, and all the gets are done when reading the fastq file.
 There are expected to be duplicate keys, since there are multiple mappings of read names to gap
 ids. An alternative would be to have each read_orient_t store an array of gap ids. However, 
 such data cannot be updated with atomics alone, so instead we just use duplicates. 
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <upc.h>
#include <upc_nb.h>

#include "utils.h"
#include "timers.h"
#include "tracing.h"
#include "dhtable.h"
#include "htable.h"

static double _load_factor = 1.0;
static int _max_depth = 0;
static long _nentries = 0;
static shared int *_entry_counts;
static shared entry_t *_entries;


void dhtable_init(long tot_reads, double load_factor) 
{
    _load_factor = load_factor;
    _nentries = tot_reads / load_factor + 3;
    //tprintf_flush("dhtable_init(%ld, %f) for %ld entries\n", tot_reads, load_factor, _nentries);
    _entry_counts = upc_all_alloc(_nentries, sizeof(int));
    _entries = upc_all_alloc(_nentries, sizeof(entry_t));
    if (_entry_counts == NULL || _entries == NULL) DIE("Could not allocate %ld entries for dhtable!\n", _nentries);
    upc_forall (long i = 0; i < _nentries; i++; i) 
        _entry_counts[i] = 0;

    UPC_TIMED_BARRIER;
}

void dhtable_free(void) 
{
    if (!MYTHREAD) {
        upc_free(_entries);
        upc_free(_entry_counts);
    }
}

// this will store duplicate entries
int dhtable_put(char *key, int value)
{
    size_t h = htable_hash(key);
    long pos = -1;
    int depth = 0;
    for (long n = _nentries / 2, offset = 0; n > 1; offset += n, n >>= 1) {
        pos = h % n + offset;
        if (!bupc_atomicI_cswap_relaxed(&_entry_counts[pos], 0, 1))
            break;
        pos = -1;
        depth++;
    }
    if (depth > _max_depth) {
		_max_depth = depth;
    }
    if (pos == -1) 
        return 0;

    upc_memput_nbi(_entries[pos].key, key, strlen(key) + 1);
    _entries[pos].value = value;
    return 1;
}

// returns the number of values and fills the list with those values
int dhtable_get(char *key, int *values, int max_values)
{
    START_TIMER(T_GET_READ_ORIENT);
    size_t h = htable_hash(key);
    long pos = -1;
    char local_key[MAX_READ_NAME_LEN + 2] = "";
    int nvals = 0;
    //char buf[200];
    for (long n = _nentries / 2, depth = 0, offset = 0; n > 1; offset += n, n >>= 1, depth++) {
        pos = h % n + offset;
        if (!_entry_counts[pos])
            break;
        upc_memget(local_key, _entries[pos].key, MAX_READ_NAME_LEN + 1);
        if (strcmp(local_key, key) == 0) {
            if (nvals == max_values)
                DIE("Too many values > %d for dhtable, key %s\n", max_values, key);
            values[nvals] = _entries[pos].value;
            nvals++;
        }
    }
	stop_timer(T_GET_READ_ORIENT);
    return nvals;
}

void dhtable_check(void) 
{
    UPC_TIMED_BARRIER;
    serial_printf("Checking read data hash...\n");
    START_TIMER(T_CHECK_DHASH);
    
    long nelems = 0;
    long max_depth = 0;
    long tot_depth = 0;
    
    const int HISTDEPTH = 30;
    int *histogram = calloc_chk(HISTDEPTH, sizeof(int));
    //char name[100];
    upc_forall (long i = 0; i < _nentries; i++; i) {
        nelems += _entry_counts[i];
        if (_entry_counts[i]) {
            //upc_memget(name, _entries[i].name, 99);
            int depth = 0;
            for (long n = _nentries / 2, offset = 0; n > 1; offset += n, n >>= 1, depth++) {
                if (i < offset + n)
                    break;
            }
            if (depth >= HISTDEPTH)
                depth = HISTDEPTH - 1;
            histogram[depth]++;
            if (depth > max_depth)
                max_depth = depth;
            tot_depth += max_depth;
        }
    }
    UPC_TIMED_BARRIER;

    long max_possible_depth = 0;
    if (!MYTHREAD) 
        for (long n = _nentries / 2, offset = 0; n > 1; offset += n, n >>= 1, max_possible_depth++);

    long all_nelems = reduce_long(nelems, UPC_ADD);
    long all_max_depth = reduce_long(max_depth, UPC_MAX);
    long all_tot_depth = reduce_long(tot_depth, UPC_ADD);

    serial_printf("\t_entries hash table stats:\n"
                  "\t  capacity:        %ld\n"
                  "\t  elements:        %ld\n"
                  "\t  expected load:   %.3f\n"
                  "\t  load:            %.3f\n"
                  "\t  max poss depth:  %ld\n"
                  "\t  max depth:       %ld\n"
                  "\t  av depth:        %.3f\n",
                  _nentries, all_nelems, _load_factor,
                  (double)all_nelems / _nentries, max_possible_depth, 
                  all_max_depth, (double)all_tot_depth / all_nelems);

    for (int i = 0; i < all_max_depth; i++) 
        serial_printf("\t  depth frac %d:    %.5f\n", 
                      i, (double)histogram[i] / (_nentries / 2 / THREADS));
    stop_timer(T_CHECK_DHASH);
    serial_printf("\tmemory used (gb):  %.3f\n"
                  "\ttime for stats:    %.3f\n",
                  (double)get_max_mem_usage_mb() / 1024, 
                  get_elapsed_time(T_CHECK_DHASH));
    UPC_TIMED_BARRIER;
}

