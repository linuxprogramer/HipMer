#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/time.h>
#include <math.h> 
#include <upc.h>
#include <string.h>
#include <assert.h>

#include "../common/common.h"
#include "meraculous.h"
#include "kmer_handling.h"
#include "packingDNAseq.h"

extern int64_t  local_allocs;

#ifndef COLL_ALLOC
/* Creates and initializes a distributed hastable - collective function */
hash_table_t* create_hash_table(int64_t  size, memory_heap_t *memory_heap, int64_t my_heap_size)
{
	hash_table_t *result;
	int64_t  i;
	int64_t  n_buckets = size;

	if (my_heap_size <= 0) my_heap_size = 1;
	
	result = (hash_table_t*) malloc_chk(sizeof(hash_table_t));
	result->size = n_buckets;
	result->table = (shared[BS] bucket_t*) upc_all_alloc(n_buckets, sizeof(bucket_t));
	if (result->table == NULL) 
        DIE("Could not allocate memory for the distributed hash table! %" PRId64 " buckets of %lu bytes\n", 
            n_buckets, (unsigned long)sizeof(bucket_t));
	
	for (i=MYTHREAD; i < n_buckets; i+=THREADS) {
		result->table[i].head = NULL;
	}
	
	memory_heap->heap_ptr = (shared[BS] shared_heap_ptr*) upc_all_alloc(THREADS, sizeof(shared_heap_ptr));
	memory_heap->heap_indices = (shared[BS] UPC_INT64_T*) upc_all_alloc(THREADS, sizeof(UPC_INT64_T));
	if (memory_heap->heap_indices == NULL || memory_heap->heap_ptr == NULL) 
        DIE("Could not allocate memory for the distributed hash table or index arrays!\n");

	memory_heap->heap_size = my_heap_size; // local variable
	memory_heap->heap_indices[MYTHREAD] = 0;
	memory_heap->heap_ptr[MYTHREAD] = (shared[] list_t*) upc_alloc( my_heap_size * sizeof(list_t));
	if (memory_heap->heap_ptr[MYTHREAD] == NULL) 
		DIE("Could not allocate memory for the distributed hash table! %" PRId64 " elements of %lu bytes\n", 
            my_heap_size, (unsigned long)sizeof(list_t));

   upc_barrier;
   
   memory_heap->cached_heap_ptrs = (shared_heap_ptr *) malloc_chk(THREADS * sizeof(shared_heap_ptr));
   memory_heap->cached_heap_indices = (shared_int64_ptr *) malloc_chk(THREADS * sizeof(shared_int64_ptr));

   for (i=0; i < THREADS; i++) {
      memory_heap->cached_heap_ptrs[i] = memory_heap->heap_ptr[i];
      memory_heap->cached_heap_indices[i] = &(memory_heap->heap_indices[i]);
   }
	
	return result;
}
#endif

#ifdef COLL_ALLOC
/* Creates and initializes a distributed hastable - collective function */
hash_table_t* create_hash_table(int64_t  size, memory_heap_t *memory_heap)
{
	hash_table_t *result;
	int64_t  i;
	int64_t  n_buckets = size;
	
#ifdef ALLOC_TIMING
	UPC_TICK_T start_timer, end_timer;
	upc_barrier;
	if (MYTHREAD == 0)
		start_timer = UPC_TICKS_NOW();
#endif

	result = (hash_table_t*) malloc_chk(sizeof(hash_table_t));
	result->size = n_buckets;
	result->table = (shared[BS] bucket_t*) upc_all_alloc(n_buckets, sizeof(bucket_t));
	
#ifdef ALLOC_TIMING
	upc_barrier;
	if (MYTHREAD == 0) {
		end_timer = UPC_TICKS_NOW();
		printf("\nupc_all_alloc()time for buckets is %f seconds\n", UPC_TICKS_TO_SECS(end_timer-start_timer);
		
	}
#endif

#ifdef ALLOC_TIMING
	if (MYTHREAD == 0)
		start_timer = UPC_TICKS_NOW();
#endif

	for (i=MYTHREAD; i < n_buckets; i+=THREADS) {
			result->table[i].head = NULL;
	}
	
#ifdef ALLOC_TIMING
	upc_barrier;
	if (MYTHREAD == 0) {
		end_timer = UPC_TICKS_NOW();
		printf("\nTime for initializing heads is %f seconds\n", UPC_TICKS_TO_SECS(end_timer-start_timer));
		
	}
#endif

#ifdef ALLOC_TIMING
	if (MYTHREAD == 0)
		start_timer = UPC_TICKS_NOW();
#endif
	shared[BS] list_t *init_ptr = (shared[BS] list_t*) upc_all_alloc(EXPANSION_FACTOR * n_buckets), sizeof(list_t));
	
#ifdef ALLOC_TIMING
	upc_barrier;
	if (MYTHREAD == 0) {
		end_timer = UPC_TICKS_NOW();
		printf("\nupc_all_alloc()time for entries is %f seconds\n", UPC_TICKS_TO_SECS(end_timer-start_timer));

	}
#endif
	
#ifdef ALLOC_TIMING
	if (MYTHREAD == 0)
		start_timer = UPC_TICKS_NOW();
#endif
	memory_heap->heap_ptr = (shared[BS] shared_heap_ptr*) upc_all_alloc(THREADS, sizeof(shared_heap_ptr));
#ifdef ALLOC_TIMING
	upc_barrier;
	if (MYTHREAD == 0) {
		end_timer = UPC_TICKS_NOW();
		printf("\nupc_all_alloc()time for heap_ptrs is %f seconds\n", UPC_TICKS_TO_SECS(end_timer-start_timer));
		
	}
#endif
	
#ifdef ALLOC_TIMING
	if (MYTHREAD == 0)
		start_timer = UPC_TICKS_NOW();
#endif

	memory_heap->heap_indices = (shared[BS] UPC_INT64_T*) upc_all_alloc(THREADS, sizeof(UPC_INT64_T));

#ifdef ALLOC_TIMING
	upc_barrier;
	if (MYTHREAD == 0) {
		end_timer = UPC_TICKS_NOW();
		printf("\nupc_all_alloc()time for heap indices is %f seconds\n", UPC_TICKS_TO_SECS(end_timer-start_timer));
		
	}
#endif

#ifdef ALLOC_TIMING
	if (MYTHREAD == 0)
		start_timer = UPC_TICKS_NOW();
#endif
	memory_heap->heap_ptr[MYTHREAD] = (shared[] list_t*)  (init_ptr + MYTHREAD);
#ifdef ALLOC_TIMING
	upc_barrier;
	if (MYTHREAD == 0) {
		end_timer = UPC_TICKS_NOW();
		printf("\nTime for initializing heap ptrs is %f seconds\n", UPC_TICKS_TO_SECS(end_timer-start_timer));
		
	}
#endif

	/* Make sure that memory is allocated before copying local object */
#ifdef ALLOC_TIMING
	if (MYTHREAD == 0)
		start_timer = UPC_TICKS_NOW();
#endif	
	memory_heap->heap_indices[MYTHREAD] = 0;
	upc_barrier;
#ifdef ALLOC_TIMING
	if (MYTHREAD == 0) {
		end_timer = UPC_TICKS_NOW();
		printf("\nTime for initializing heap indices is %f seconds\n", UPC_TICKS_TO_SECS(end_timer-start_timer));
		
	}
#endif
	
	return result;
}
#endif

/* Frees the distributed hashtable */
/*void free_table(hash_table_t *hashtable, memory_heap_t *memory_heap)
{
    int64_t  i;	
    if (hashtable==NULL) return;
	
    upc_forall(i=0; i < hashtable->size; i++; i) {
#ifndef STORE_OPT
		upc_lock_free(hashtable->table[i].bucket_lock);
#endif
    }
	
	upc_barrier;
	upc_free(hashtable->table);
	//upc_free(memory_heap->heap_ptr);
	upc_free(memory_heap->heap_indices);
	
    free(hashtable);
	upc_barrier;
}*/

/* Computes the hash value of null terminated sequence */
/*
int64_t  hash(int64_t  hashtable_size, char *kmer)
{
    unsigned long hashval;
    hashval = 5381;
    for(; *kmer != '\0'; kmer++) hashval = (*kmer) +  (hashval << 5) + hashval;
      
    return hashval % hashtable_size;
}

uint64_t  hashseq(int64_t  hashtable_size, char *seq, int size)
{
	unsigned long hashval;
	hashval = 5381;
	for(int i = 0; i < size; i++) {
		hashval = seq[i] +  (hashval << 5) + hashval;
	}
	
	return hashval % hashtable_size;
}
*/

uint64_t hashkmer(int64_t  hashtable_size, char *seq)
{
	return hashkey(hashtable_size, seq, KMER_LENGTH);
}
   
#ifdef LSH

/* Function that loads oracle table from file to node's shared memory */
int loadOracleTable(int64_t nBucketsIn, char *oracleFile ,shared[1] oracle_t **result)
{
   int64_t nBuckets = (nBucketsIn / cores_per_node) * cores_per_node + cores_per_node ;
   int64_t nodes = (int64_t) ceil((double) THREADS / (double)cores_per_node);
   int64_t oracleTablesSize = nodes * nBuckets;
   FILE *fd = fopen_chk( oracleFile , "r");
   
   shared[1] oracle_t *tmpTable;

   tmpTable = (shared[1] oracle_t*) upc_all_alloc(oracleTablesSize , sizeof(oracle_t));
   
   int64_t myInId = MYTHREAD % cores_per_node;
   int64_t offset = nBucketsIn / cores_per_node * myInId;
   fseek(fd, offset * sizeof(oracle_t), SEEK_SET);
   //fseek(fd, offset, SEEK_SET);
   
   int64_t my_lines;
   if (myInId == (cores_per_node-1)) {
      my_lines = (nBucketsIn/cores_per_node + nBucketsIn%cores_per_node);
   } else {
      my_lines = nBucketsIn/cores_per_node;
   }
   
   int64_t cur_line;
   int cur_node;
   int64_t loc;
   int64_t nid = MYTHREAD / cores_per_node;
   int64_t chunk;
   int64_t offset2;
   oracle_t *input = (oracle_t*) malloc_chk(my_lines*sizeof(oracle_t));
   fread(input, sizeof(oracle_t), my_lines, fd);
   
   for (cur_line = 0; cur_line < my_lines; cur_line++) {
      chunk = (offset+cur_line) / cores_per_node;
      loc = nid * cores_per_node + chunk * THREADS + ((offset+cur_line) % cores_per_node);
      tmpTable[loc] = (oracle_t) input[cur_line];
   }
   
   (*result) = tmpTable;
   
   return 1;
}

int64_t LSHkmer(int64_t hashtable_size, char *seq)
{
   
	int64_t original_hashval = hashkmer(hashtable_size, seq);
   int64_t chunk = original_hashval / cores_per_node;
   int64_t offset = original_hashval % cores_per_node;
   int64_t nid = MYTHREAD / cores_per_node;
   int64_t nodes = (int64_t) ceil((double) THREADS / (double)cores_per_node);
   int64_t loc = nid * cores_per_node + chunk * THREADS + offset;
   oracle_t node_id = (oracle_t) oracle_table[loc];
   
   int64_t repetitive_bucket = original_hashval / THREADS;
   int64_t newVal = repetitive_bucket * THREADS + node_id * cores_per_node + offset;
   
   return newVal;
}
   
#endif
   
/*   int64_t LSHkmer(int64_t hashtable_size, char *seq)
   {
      int shingleLen = 14;
      int nShingles = KMER_LENGTH - shingleLen + 1;
      int i;
      
      int64_t min = hashseq(hashtable_size, seq, shingleLen);
      int64_t newVal;
      
      for (i = 1; i < nShingles; i++) {
         newVal = hashseq(hashtable_size, seq+i, shingleLen);
         min = (newVal < min) ? newVal : min;
      }
      
      return min;
   } */
 
/*int64_t  hash(int hashtable_size, char *kmer)
{
    int64_t  hashval, i;
		
	for (hashval=KMER_LENGTH, i=0; i<KMER_LENGTH; ++i)
		hashval = (hashval<<4)^(hashval>>28)^kmer[i];
	return ((hashval % 1610612741) % hashtable_size);
}*/

/* Use this lookup function when no writes take place in the distributed hashtable */
shared[] list_t* lookup_kmer_and_copy(hash_table_t *hashtable, const unsigned char *kmer, list_t *cached_copy)
{	
	// TODO: Fix hash functions to compute on packed kmers and avoid conversions
#ifdef LSH
   int64_t  hashval = LSHkmer(hashtable->size, (char*) kmer);
#else
   int64_t  hashval = hashkmer(hashtable->size, (char*) kmer);
#endif
	unsigned char packed_key[KMER_PACKED_LENGTH];
	unsigned char remote_packed_key[KMER_PACKED_LENGTH];

	packSequence(kmer, packed_key, KMER_LENGTH);
	shared[] list_t *result;
	bucket_t local_buc;
	
	local_buc = hashtable->table[hashval];
	result = local_buc.head;
	

	for (; result != NULL;) {
		//upc_memget(remote_packed_key, result->packed_key, KMER_PACKED_LENGTH * sizeof(char));
		//upc_memget(local_res, result, sizeof(list_t));
		assert( upc_threadof( result ) == hashval % THREADS );
#ifdef MERACULOUS
		loop_until( *cached_copy = *result , IS_VALID_UPC_PTR(cached_copy->my_contig) );
#endif
#ifdef MERDEPTH
      *cached_copy = *result;
#endif
#ifdef CEA
      *cached_copy = *result;
#endif
		//if (local_res.structid != 1988) {
		//	printf("FATAL error in lookup, struct id is %d!!!\n", local_res.structid);
		//}
		if (comparePackedSeq(packed_key, cached_copy->packed_key, KMER_PACKED_LENGTH) == 0){
			return result;
		}
		result = cached_copy->next;
	}
	if (VERBOSE > 2) {
		unsigned char tmp[KMER_LENGTH+1];
		tmp[KMER_LENGTH] = '\0';
		unpackSequence(packed_key, tmp, KMER_LENGTH);
		LOG("Thread %d: lookup_kmer(): did not find %s (%s)\n", MYTHREAD, tmp, kmer);
	}
	return NULL;
}

shared[] list_t* lookup_kmer(hash_table_t *hashtable, const unsigned char *kmer) {
	list_t cached_copy;
	return lookup_kmer_and_copy(hashtable, kmer, &cached_copy);
}

/* find the entry for this kmer or reverse complement */
shared[] list_t *lookup_least_kmer_and_copy(hash_table_t *dist_hashtable, const char *next_kmer, list_t *cached_copy, int *is_least) {
	char auxiliary_kmer[KMER_LENGTH+1];
	char *kmer_to_search;
	auxiliary_kmer[KMER_LENGTH] = '\0';
	/* Search for the canonical kmer */
	kmer_to_search = getLeastKmer(next_kmer, auxiliary_kmer);
	*is_least = (kmer_to_search == next_kmer) ? 1 : 0;
	shared[] list_t *res = lookup_kmer_and_copy(dist_hashtable, (unsigned char*) kmer_to_search, cached_copy);
	if (VERBOSE > 2) LOG("Thread %d: lookup_least_kmer2(%s) is_least: %d res: %s\n", MYTHREAD, next_kmer, *is_least, res == NULL ? " not found" : "found");
	return res;
}

shared[] list_t *lookup_least_kmer(hash_table_t *dist_hashtable, const char *next_kmer, int *was_least) {
	list_t cached_copy;
	return lookup_least_kmer_and_copy(dist_hashtable, next_kmer, &cached_copy, was_least);
}

void set_ext_of_kmer(list_t *lookup_res, int is_least, char *new_seed_le, char *new_seed_re) {
	if (is_least) {
#ifdef MERACULOUS
		convertPackedCodeToExtension(lookup_res->packed_extensions,new_seed_le,new_seed_re);
#endif
	} else {
#ifdef MERACULOUS
		convertPackedCodeToExtension(lookup_res->packed_extensions,new_seed_re,new_seed_le);
#endif
		*new_seed_re = reverseComplementBaseExt(*new_seed_re);
		*new_seed_le = reverseComplementBaseExt(*new_seed_le);
	}
	if (VERBOSE > 2) LOG("Thread %d: set_ext_of_kmer: is_least: %d l:%c r:%c\n", MYTHREAD, is_least, *new_seed_le, *new_seed_re);
	
}

/* find the entry for this kmer or reverse complement, set the left and right extensions */
shared[] list_t *lookup_and_get_ext_of_kmer(hash_table_t *dist_hashtable, const char *next_kmer, char *new_seed_le, char *new_seed_re)
{
	int is_least;
	shared[] list_t *lookup_res;
	list_t copy;
	*new_seed_le = '\0';
	*new_seed_re = '\0';
	
	/* Search for the canonical kmer */
	lookup_res = lookup_least_kmer_and_copy(dist_hashtable, next_kmer, &copy, &is_least);
	
	if (lookup_res == NULL) {
		return lookup_res;
	}
	
	/* Find extensions of the new kmer found in the hashtable */
	set_ext_of_kmer(&copy, is_least, new_seed_le, new_seed_re);
	
	return lookup_res;
}

/* get the contig_ptr "box" associated with a kmer -- must be USED already to work */
shared[] contig_ptr_box_list_t *get_contig_box_of_kmer(shared[] list_t *lookup_res, int follow_list) {
	shared[] contig_ptr_box_list_t *box_ptr = NULL;
	kmer_and_ext_t kmer_and_ext;
	
	if (lookup_res == NULL)
		printf("FATAL ERROR: K-mer should not be NULL here (right walk)\n");
	assert(lookup_res != NULL);
#ifdef MERACULOUS
	assert(lookup_res->used_flag == USED);
	loop_until ( box_ptr = lookup_res->my_contig, box_ptr != NULL && IS_VALID_UPC_PTR(box_ptr)  );
#endif
	// follow box to tail
	while(follow_list && box_ptr->next != NULL) {
		box_ptr = box_ptr->next;
	}

	return box_ptr;
}

/* get the contig associated with a kmer -- must be USED already to work */
shared[] contig_t *get_contig_of_kmer(shared[] list_t *lookup_res, int follow_list) {
	shared[] contig_ptr_box_list_t *box_ptr;
	kmer_and_ext_t kmer_and_ext;
	shared[] contig_t *contig = NULL;
	
	assert(lookup_res != NULL);
#ifdef MERACULOUS
	assert(lookup_res->used_flag == USED);
#endif

	loop_until( 
		box_ptr = get_contig_box_of_kmer(lookup_res, follow_list); contig = box_ptr->contig,
		contig != NULL && IS_VALID_UPC_PTR(contig)
		);

	// can not make this assertion in parallel as the list can grow
	// assert( follow_list == 0 || box_ptr->next == NULL );
	
	return contig;
}

/* Use this lookup function when writes take place in the distributed hashtable - i.e. during creation */
/*shared[] list_t* sync_lookup_kmer(hash_table_t *hashtable, unsigned char *kmer)
{	
	// TODO: Fix hash functions to compute on packed kmers and avoid conversions
    int64_t  hashval = hash(hashtable->size, (char*) kmer);
	unsigned char packed_key[KMER_PACKED_LENGTH];
	packSequence(kmer, packed_key, KMER_LENGTH);
	shared[] list_t *result;
	list_t remote_element;
	
	upc_lock(hashtable->table[hashval].bucket_lock);
	for (result = hashtable->table[hashval].head; result != NULL; result = result->next) {
		upc_memget(&remote_element, result, sizeof(list_t));
		if (comparePackedSeq(packed_key, remote_element.packed_key, KMER_PACKED_LENGTH) == 0){ 
			upc_unlock(hashtable->table[hashval].bucket_lock);
			return result;
		}
	}
	upc_unlock(hashtable->table[hashval].bucket_lock);
	return NULL;
}*/

#ifndef AGGREGATE_STORES
int add_kmer(hash_table_t *hashtable, unsigned char *kmer, unsigned char *extensions, memory_heap_t *memory_heap, int64_t  *ptrs_to_stack, int64_t  *ptrs_to_stack_limits, int CHUNK_SIZE)
{
	list_t new_entry;
	shared[] list_t *insertion_point;

#ifdef LSH
	int64_t  hashval = LSHkmer(hashtable->size, (char*) kmer);
#else
   int64_t  hashval = hashkmer(hashtable->size, (char*) kmer);
#endif
	int64_t  pos;
	int64_t  new_init_pos, new_limit;
	int remote_thread = upc_threadof(&(hashtable->table[hashval]));
	
	if (ptrs_to_stack[remote_thread] < ptrs_to_stack_limits[remote_thread]) {
		pos = ptrs_to_stack[remote_thread];
		ptrs_to_stack[remote_thread]++;
	} else {
	
		new_init_pos = UPC_ATOMIC_FADD_I64(&(memory_heap->heap_indices[remote_thread]), CHUNK_SIZE);
		
		new_limit = new_init_pos + CHUNK_SIZE;
		if (new_limit > EXPANSION_FACTOR * (hashtable->size / THREADS))
			new_limit = EXPANSION_FACTOR * (hashtable->size / THREADS);
		if (new_limit <= new_init_pos)
			pos = EXPANSION_FACTOR *(hashtable->size / THREADS);
		else
			pos = new_init_pos;
		new_init_pos++;
		ptrs_to_stack[remote_thread] = new_init_pos;
		ptrs_to_stack_limits[remote_thread] = new_limit;
	}
	
	if (pos < EXPANSION_FACTOR * (hashtable->size / THREADS)) {
		// Case 1: There is more space in the remote heap
		insertion_point = (shared[] list_t*) ((memory_heap->heap_ptr[remote_thread]) + pos);
	} else {
		// Case 2; There is NOT more space in the remote heap, allocate locally space
		insertion_point = (shared[] list_t*) upc_alloc(sizeof(list_t));
#ifdef MERACULOUS
		insertion_point->my_contig  = (shared[] contig_ptr_box_list_t*) upc_alloc(sizeof(contig_ptr_box_list_t));
		insertion_point->my_contig->contig = NULL;
#endif
		//local_allocs++;
	}
	
	/* Pack and insert kmer info into list - use locks to avoid race conditions */
	packSequence(kmer, new_entry.packed_key, KMER_LENGTH);
#ifdef MERACULOUS
	new_entry.packed_extensions = convertExtensionsToPackedCode(extensions);
	new_entry.used_flag = UNUSED;
#endif
	new_entry.next = NULL;

#ifdef SYNC_PROTOCOL
#ifdef MERACULOUS
	new_entry.my_contig = NULL;
#endif
#endif
	upc_memput(insertion_point, &new_entry, sizeof(list_t));
	//printf("Thread %d returned from upc_memput() and hashval is %lld \n", MYTHREAD, hashval);

	/* Use locks to update the head pointer in the remote bucket */
#ifndef STORE_OPT
	upc_lock(hashtable->table[hashval].bucket_lock);
#endif
	insertion_point->next = hashtable->table[hashval].head;
	hashtable->table[hashval].head = insertion_point;
#ifndef STORE_OPT
	upc_unlock(hashtable->table[hashval].bucket_lock);
#endif
	//printf("Thread %d completed add_kmer for hashval %lld \n", MYTHREAD, hashval);

	
	return 0;
}
#endif

#ifdef AGGREGATE_STORES
int add_kmer(hash_table_t *hashtable, unsigned char *kmer, unsigned char *extensions, memory_heap_t *memory_heap, int64_t  *local_buffs_hashvals, list_t *local_buffs, int64_t  *local_index, int CHUNK_SIZE)
{
   const double EXPANSION_FACTOR = 1.2;

	list_t new_entry;
	shared[] list_t *insertion_point;
#ifdef LSH
   int64_t  hashval = LSHkmer(hashtable->size, (char*) kmer);
#else
   int64_t  hashval = hashkmer(hashtable->size, (char*) kmer);
#endif	
   int64_t  j, excess;
	int64_t  store_pos;
	int remote_thread = upc_threadof(&(hashtable->table[hashval]));
	int64_t  cur_hashval;
	
	/* Pack and insert kmer info into local_list_t buffer */
	packSequence(kmer, new_entry.packed_key, KMER_LENGTH);
	new_entry.packed_extensions = convertExtensionsToPackedCode(extensions);
#ifdef MERACULOUS
	new_entry.used_flag = UNUSED;
#endif
	new_entry.next = NULL;
#ifdef SYNC_PROTOCOL
	new_entry.my_contig = NULL;
#endif
	
	memcpy( &(local_buffs[local_index[remote_thread] + remote_thread * CHUNK_SIZE]), &new_entry, sizeof(list_t));
	local_buffs_hashvals[local_index[remote_thread] + remote_thread * CHUNK_SIZE] = hashval;
	
	if (local_index[remote_thread] == CHUNK_SIZE - 1) {        
		/* Now we should do a batch upc_memput */
		store_pos = UPC_ATOMIC_FADD_I64(&(memory_heap->heap_indices[remote_thread]), CHUNK_SIZE);

		if (store_pos >= (EXPANSION_FACTOR * ((hashtable->size + hashtable->size % THREADS) / THREADS))) {
			store_pos = (EXPANSION_FACTOR * ((hashtable->size + hashtable->size % THREADS) / THREADS));
		}
		
		/* Check if there is sufficient space left in the remote heap  */
		if ((store_pos + CHUNK_SIZE) > (EXPANSION_FACTOR * ((hashtable->size + hashtable->size % THREADS) / THREADS))) {
			excess = (store_pos + CHUNK_SIZE) - (EXPANSION_FACTOR * ((hashtable->size + hashtable->size % THREADS) / THREADS));
			upc_memput( (shared list_t*)  ((memory_heap->heap_ptr[remote_thread]) + store_pos)  , &(local_buffs[remote_thread * CHUNK_SIZE]), (CHUNK_SIZE-excess) * sizeof(list_t));
			for (j = 0; j < CHUNK_SIZE - excess; j++) {
				/* Use locks to update the head pointers in the remote buckets */
				cur_hashval = local_buffs_hashvals[j + remote_thread * CHUNK_SIZE];
				upc_lock(hashtable->table[cur_hashval].bucket_lock);
				((memory_heap->heap_ptr[remote_thread]) + store_pos + j)->next = hashtable->table[cur_hashval].head;
				hashtable->table[cur_hashval].head = (shared list_t*) ((memory_heap->heap_ptr[remote_thread]) + store_pos + j);
				upc_unlock(hashtable->table[cur_hashval].bucket_lock);
			}
			insertion_point = (shared[] list_t*) upc_alloc(excess * sizeof(list_t));
			upc_memput( insertion_point , &(local_buffs[remote_thread * CHUNK_SIZE + (CHUNK_SIZE-excess)]), excess * sizeof(list_t));
			for (j = 0; j < excess; j++) {
				/* Use locks to update the head pointers in the remote buckets */
				cur_hashval = local_buffs_hashvals[CHUNK_SIZE - excess + j + remote_thread * CHUNK_SIZE];
				upc_lock(hashtable->table[cur_hashval].bucket_lock);
				(insertion_point + j)->next = hashtable->table[cur_hashval].head;
				hashtable->table[cur_hashval].head = (shared list_t*) (insertion_point + j);
				upc_unlock(hashtable->table[cur_hashval].bucket_lock);
			}
			
		} else {
			upc_memput( (shared list_t*)  ((memory_heap->heap_ptr[remote_thread]) + store_pos)  , &(local_buffs[remote_thread * CHUNK_SIZE]), CHUNK_SIZE * sizeof(list_t));
			for (j = 0; j < CHUNK_SIZE; j++) {
				/* Use locks to update the head pointers in the remote buckets */
				cur_hashval = local_buffs_hashvals[j + remote_thread * CHUNK_SIZE];
				upc_lock(hashtable->table[cur_hashval].bucket_lock);
				((memory_heap->heap_ptr[remote_thread]) + store_pos + j)->next = hashtable->table[cur_hashval].head;
				hashtable->table[cur_hashval].head = (shared list_t*) ((memory_heap->heap_ptr[remote_thread]) + store_pos + j);
				upc_unlock(hashtable->table[cur_hashval].bucket_lock);
			}        
		}
		
		local_index[remote_thread] = 0;
	} else {
		local_index[remote_thread]++;
	}
	
    return 0;
}
#endif

/*
int delete_kmer(hash_table_t *hashtable, unsigned char *kmer)
{
	int64_t i;
	shared[] list_t *list, *prev;
	int64_t  hashval = hash(hashtable->size, (char*) kmer);
	unsigned char packed_key[KMER_PACKED_LENGTH];
	packSequence(kmer, packed_key, KMER_LENGTH);
	list_t to_be_deleted;
	
	upc_lock(hashtable->table[hashval].bucket_lock);
	for(prev=NULL, list=hashtable->table[hashval].head; list != NULL; prev = list, list = list->next) {
		upc_memget(&to_be_deleted, list, sizeof(list_t));
		if (comparePackedSeq(packed_key, to_be_deleted.packed_key, KMER_PACKED_LENGTH) == 0) break;
	}
	
	if (list==NULL) return 1; 
	
	if (prev==NULL) hashtable->table[hashval].head = list->next;
	else prev->next = list->next; 
	
	//upc_free(list);	
	
	upc_unlock(hashtable->table[hashval].bucket_lock);
	
	return 0;
}*/

#endif // KMER_HASH_H

