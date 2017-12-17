#ifndef __KMER_HASH2_H
#define __KMER_HASH2_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <upc.h>
#include <string.h>
#include "packingDNAseq2.h"
#include "merAligner.h"
#include "../common/common.h"

extern int64_t  local_allocs;

/* Creates and initializes a distributed hastable - collective function */
hash_table_t* create_hash_table(int64_t size, memory_heap_t *memory_heap, double EXPANSION_FACTOR)
{
	hash_table_t *result;
	int64_t  i;
	int64_t  n_buckets = size;
	int64_t  expanded_size;
	
	result = (hash_table_t*) malloc_chk(sizeof(hash_table_t));
	result->size = n_buckets;
	result->table = (shared[BS] bucket_t*) upc_all_alloc(n_buckets, sizeof(bucket_t));
	
	for (i=MYTHREAD; i < n_buckets; i+=THREADS) {
		assert( upc_threadof( &(result->table[i]) ) == MYTHREAD );
		result->table[i].head = NULL;
#ifndef STORE_OPT
		result->table[i].bucket_lock = upc_global_lock_alloc();
#endif
	}
	
	memory_heap->heap_ptr = (shared[BS] shared_heap_ptr*) upc_all_alloc(THREADS, sizeof(shared_heap_ptr));
	memory_heap->heap_indices = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
	memory_heap->contigInfoArrays = (shared[BS] contigInfoArrayPtr_t *) upc_all_alloc(THREADS, sizeof(contigInfoArrayPtr_t));
        memory_heap->heap_ptr_sizes = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
        expanded_size = ( (int64_t) (ceil(EXPANSION_FACTOR * (size + size % THREADS)))) / THREADS * sizeof(list_t);
    if (memory_heap->heap_ptr == NULL || memory_heap->heap_indices == NULL || memory_heap->contigInfoArrays == NULL || memory_heap->heap_ptr_sizes == NULL) 
        DIE("Failed to allocate a heap and indicies for hashtable\n");
	memory_heap->heap_ptr[MYTHREAD] = (shared[] list_t*) upc_alloc(expanded_size);
    if (memory_heap->heap_ptr[MYTHREAD] == NULL) 
        DIE("Failed to allocate a heap for hashtable %ld\n", expanded_size); 
	memory_heap->heap_indices[MYTHREAD] = 0;
        memory_heap->contigInfoArrays[MYTHREAD] = NULL;
        memory_heap->heap_ptr_sizes[MYTHREAD] = expanded_size;
	upc_barrier;
	
	return result;
}

/* Computes the hash value of a k-mer key */
/*
int64_t  hash(int64_t  hashtable_size, char *kmer)
{
    unsigned long hashval;
    hashval = 5381;
    for(; *kmer != '\0'; kmer++) hashval = (*kmer) +  (hashval << 5) + hashval;
      
    return hashval % hashtable_size;
}
*/

shared[] list_t* lookup_kmer_in_bucket(shared[] bucket_t *bucket, const unsigned char *packed_key) {
	shared[] list_t *result = bucket->head;
	list_t local_res;
	while(result != NULL) {
		local_res = *result;
		if (comparePackedSeq(packed_key, local_res.packed_key, KMER_PACKED_LENGTH) == 0) {
			return result;
		}
		result = local_res.next;
	}
	return NULL;
}
		
		
/* Use this lookup function when no writes take place in the distributed hashtable */
shared[] list_t* lookup_kmer(hash_table_t *hashtable, const unsigned char *kmer)
{	
	int64_t  hashval = hashstr(hashtable->size, (char*) kmer);
	unsigned char packed_key[KMER_PACKED_LENGTH];
	unsigned char remote_packed_key[KMER_PACKED_LENGTH];

	packSequence(kmer, packed_key, KMER_LENGTH);
	shared[] bucket_t *bucket = (shared[] bucket_t *) &(hashtable->table[hashval]);
	return lookup_kmer_in_bucket(bucket, packed_key);
}

// returns the entry from the bucket that was either already present or this new one
shared[] list_t *add_entry_to_bucket(shared[] bucket_t *bucket, shared[] list_t *new_entry) {
#ifdef STORE_OPT
	assert( upc_threadof(new_entry) == MYTHREAD );
	assert( upc_threadof(bucket) == MYTHREAD );
	new_entry->next = bucket->head;
	bucket->head = new_entry;
	return new_entry;
#else
	upc_lock(bucket->bucket_lock);
	list_t localEntry = *new_entry;
	bucket_t localBucket = *bucket;
	shared[] list_t *oldHead = localBucket.head;
	shared[] list_t *cur_ptr = lookup_kmer_in_bucket(bucket, localEntry.packed_key);
	
	if (cur_ptr == NULL) {
		// add to the list
		assert(bucket->head == oldHead);
		assert(new_entry->next == NULL);
		new_entry->next = oldHead;
		bucket->head = new_entry;
		cur_ptr = new_entry;
        	upc_fence;
	}
        upc_unlock(bucket->bucket_lock);
	return cur_ptr;
#endif
}

#endif

