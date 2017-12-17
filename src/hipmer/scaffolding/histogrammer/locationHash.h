#ifndef __LOCATION_HASH_H
#define __LOCATION_HASH_H

#include "../../common/common.h"

#define BS 1
#define EXPANSION_FACTOR 2

/* Location data type */
typedef struct location_t location_t;
struct location_t{
   int subjectID;
   int loc1;
   int loc2;
   int count;
   shared[] location_t *next;						// Pointer to next entry in the same bucket
};

/* Bucket data structure */
typedef struct bucket_t bucket_t;
struct bucket_t{
	upc_lock_t *bucket_lock;						// Lock to avoid race conditions during updates
	shared[] location_t *head;						// Pointer to the first entry of the hashtable
};

/* Hash-table data structure */
typedef struct hash_table_t hash_table_t;
struct hash_table_t {
	int64_t size;										// Size of the hashtable
	shared[BS] bucket_t *table;					// Entries of the hashtable are pointers to buckets
};

typedef shared[] location_t* shared_heap_ptr;

/* Memory heap data structure */
typedef struct memory_heap_t memory_heap_t;
struct memory_heap_t {
	shared[BS] shared_heap_ptr *heap_ptr;					// Pointers to shared memory heaps
	shared[BS] int64_t *heap_indices;						// Indices of remote heaps
        int64_t heap_size;
};

typedef struct contigScaffoldMap_t contigScaffoldMap_t;
struct contigScaffoldMap_t {
   int cStrand;
   int scaffID;
   int sStart;
   int sEnd;
};

/* Creates and initializes a distributed hash table - collective function */
hash_table_t* create_hash_table(int64_t size, memory_heap_t *memory_heap)
{
	if (MYTHREAD == 0) fprintf(stderr, "location_t hashtable size: %lld\n", (long long) size);
	hash_table_t *result;
	size = size + THREADS * 16;
	int64_t i;
	int64_t n_buckets = size;
	
	result = (hash_table_t*) malloc_chk(sizeof(hash_table_t));
	result->size = n_buckets;
	result->table = (shared[BS] bucket_t*) upc_all_alloc(n_buckets, sizeof(bucket_t));
	
	for (i=MYTHREAD; i < n_buckets; i+=THREADS) {
		result->table[i].head = NULL;
		result->table[i].bucket_lock = upc_global_lock_alloc();
	}
	
	memory_heap->heap_ptr = (shared[BS] shared_heap_ptr*) upc_all_alloc(THREADS, sizeof(shared_heap_ptr));
	memory_heap->heap_indices = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
	memory_heap->heap_size = EXPANSION_FACTOR * (size + THREADS) / THREADS + 8192;
	memory_heap->heap_ptr[MYTHREAD] = (shared[] location_t*) upc_alloc( memory_heap->heap_size * sizeof(location_t) );
	memory_heap->heap_indices[MYTHREAD] = 0;
	upc_barrier;
	
	return result;
}

/* Computes the hash value of a location entry */
/*
int64_t hash(int64_t hashtable_size, char *location)
{
   unsigned long hashval;
   hashval = 5381;
   for(; *location != '\0'; location++) hashval = (*location) +  (hashval << 5) + hashval;
   
   return hashval % hashtable_size;
}
*/

/* Look up a location in the hash table. If found return the count and increment that. Otherwise insert that location entry. This functionos thread safe, thus it holds the lock while operating on a bucket */
int lookupAndIncrementLocation(int subjectID, int loc1, int loc2, hash_table_t *hashtable, memory_heap_t *memory_heap)
{
   char locationKey[30];
   int64_t hashval;
   int count;
   int foundLocation = 0;
   shared[] location_t *result;
   location_t local_res;
   int64_t store_pos;
   int remote_thread;
   int emptyBucket = 0;
   location_t new_entry;
   shared[] location_t *heapPosition;

   
   sprintf(locationKey, "%d.%d.%d", subjectID, loc1, loc2);
   hashval = hashstr(hashtable->size, (char*) locationKey);
   
   upc_lock(hashtable->table[hashval].bucket_lock);
   
   result = hashtable->table[hashval].head;
   emptyBucket = (result == NULL) ? 1 : 0;
   for (; (result != NULL) && (foundLocation == 0);) {
      local_res = *result;
      if ((local_res.subjectID == subjectID) && (local_res.loc1 == loc1) && (local_res.loc2 == loc2)) {
         result->count++;
         count = local_res.count;
         foundLocation = 1;
      } else {
         result = local_res.next;
      }
   }
   
   if (foundLocation == 0) {
      count = 0;
      remote_thread = upc_threadof(&(hashtable->table[hashval]));
      store_pos = bupc_atomicI64_fetchadd_strict(&(memory_heap->heap_indices[remote_thread]), 1);
      if ( store_pos >= memory_heap->heap_size ) { 
          fprintf(stderr, "Thread %d: overrun on location_t heap %d: %lld over %lld\n", 
                  MYTHREAD, remote_thread, (long long) store_pos, (long long) memory_heap->heap_size); 
          upc_global_exit(1); 
      }
      new_entry.subjectID = subjectID;
      new_entry.loc1 = loc1;
      new_entry.loc2 = loc2;
      new_entry.count = 1;
      if (emptyBucket) {
         new_entry.next = NULL;
      } else {
         new_entry.next = hashtable->table[hashval].head;
      }
      heapPosition = (shared[] location_t*)  ((memory_heap->heap_ptr[remote_thread]) + store_pos);
      upc_memput(heapPosition , &new_entry, sizeof(location_t));
      hashtable->table[hashval].head = heapPosition;
   }
   
   upc_unlock(hashtable->table[hashval].bucket_lock);
   
   return count;
}

#endif
