#ifndef BUILD_UFX_HASH_H
#define BUILD_UFX_HASH_H

#define LINE_SIZE (KMER_LENGTH+4)
#include <stdio.h>
#include <sys/stat.h>
#include <assert.h>

#include "meraculous.h"
#include "packingDNAseq.h"
#include "../contigs/kmer_hash.h"

hash_table_t* buildUFXhash(int64_t size, FILE *fd, memory_heap_t *memory_heap_res, int64_t myShare, int64_t dsize, int64_t dmin, int CHUNK_SIZE, int load_factor)
{
	hash_table_t *dist_hashtable;
	memory_heap_t memory_heap;
	shared list_t *lookup_res;
	int64_t chars_read, cur_chars_read , ptr, retval, i;
	int64_t buffer_size, offset;
	int64_t num_chars;
	shared[1] int64_t *heap_sizes;
	int64_t *my_heap_sizes, my_heap_size, my_ufx_lines, idx;
	int *ufx_remote_thread, remote_thread;
	int64_t  hashval;

	UPC_TICK_T start_read, end_read, start_storing, end_storing, start_setup, end_setup, start_calculation, end_calculation;

	heap_sizes = (shared[1] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
	heap_sizes[MYTHREAD] = 0;
	
	//my_ufx_lines = (size / THREADS + size % THREADS);
	//if (MYTHREAD == THREADS-1) {
	//	my_ufx_lines = (size / THREADS + size % THREADS);
	//} else {
	//	my_ufx_lines = size / THREADS;
	//}
   
   my_ufx_lines = myShare;
	
	ufx_remote_thread = (int*) malloc_chk(my_ufx_lines * sizeof(int));
	
	if (VERBOSE > 1) LOG("Thread %d: Preparing to read %lld UFX lines\n", MYTHREAD, (long long) my_ufx_lines);
	
#ifdef DEBUG
	int64_t kmer_count = 0;
#endif
	
	/* Initialize lookup-table --> necessary for packing routines */
	init_LookupTable();
   char **kmersarr;
   int *counts;
   char *lefts;
   char *rights;
	
   start_read = UPC_TICKS_NOW();

   int64_t kmers_read = UFXRead(dsize, &kmersarr, &counts, &lefts, &rights, my_ufx_lines, dmin, 0, MYTHREAD);

   if (VERBOSE > 1) LOG("Thread %d: Loaded %lld kmers\n", MYTHREAD, (long long) kmers_read);
   
   //FILE *aux = fopen_chk("conerted.txt", "w");
   //char curkmer[KMER_LENGTH+1];
   //curkmer[KMER_LENGTH] = '\0';
   //for (i=0; i< kmers_read; i++) {
   //   fprintf( aux ,"%s\t%c%c\n", kmersarr[i], lefts[i], rights[i]);

   //}
	
	upc_barrier;
	if (MYTHREAD == 0)
		printf("Threads done with I/O\n");
	
	end_read = UPC_TICKS_NOW();
   fileIOTime = UPC_TICKS_TO_SECS(end_read-start_read);

	my_heap_sizes = (int64_t*) calloc_chk(THREADS, sizeof(int64_t));
   
   start_calculation= UPC_TICKS_NOW();
	/* calculate the exact number of UFX entries for each thread */
	/* previously estimated as EXPANSION_FACTOR * (size + size % THREADS) / THREADS */
	idx = 0;
   char rc_kmer[KMER_LENGTH+1];
   rc_kmer[KMER_LENGTH] = '\0';
   int is_least;
   
   for (ptr= 0; ptr < kmers_read; ptr++) {
      reverseComplementKmer(kmersarr[ptr], rc_kmer);
      is_least = 0;
      is_least = (strcmp(kmersarr[ptr], rc_kmer) > 0) ? 0  : 1 ;

#ifdef MERACULOUS
      if (is_least && (lefts[ptr] != 'F') && (lefts[ptr] != 'X') && (rights[ptr] != 'F') && (rights[ptr] != 'X'))
#endif
         
#ifdef MERDEPTH
      if (is_least)
#endif
         
#ifdef CEA
      if (is_least)
#endif
      {
#ifdef LSH
         hashval = LSHkmer(load_factor * size, (char*) kmersarr[ptr]);
#else
         hashval = hashkmer(load_factor * size, (char*) kmersarr[ptr]);
#endif
         remote_thread = hashval % (THREADS*BS);
         my_heap_sizes[remote_thread]++;
         ufx_remote_thread[idx] = remote_thread;
      } else {
         ufx_remote_thread[idx] = -1;
      }
      idx++;
   }
   
	assert(idx == my_ufx_lines);
	
	/* now atomically add all my_sizes to global size heap */
	/* TODO replace with all reduce? */
	for(i = MYTHREAD; i < THREADS+MYTHREAD; i++) {
		if (VERBOSE > 1) LOG("Thread %d: sending %ld to thread %d\n", MYTHREAD, my_heap_sizes[i%THREADS], (int) (i%THREADS));
		assert( upc_threadof( heap_sizes + (i%THREADS) ) == (i%THREADS) );
		UPC_ATOMIC_FADD_I64(heap_sizes + (i%THREADS), my_heap_sizes[i%THREADS]);
	}
	free(my_heap_sizes);
	upc_barrier;
	assert(upc_threadof( heap_sizes + MYTHREAD) == MYTHREAD);
	my_heap_size = heap_sizes[MYTHREAD];
	upc_all_free(heap_sizes); // should clean up after last thread gets here
	if (VERBOSE > 1) LOG("Thread %d: allocating memory for %ld kmers\n",MYTHREAD, my_heap_size);
   
   end_calculation= UPC_TICKS_NOW();
   cardCalcTime = UPC_TICKS_TO_SECS(end_calculation-start_calculation);
   
	start_setup = UPC_TICKS_NOW();

	/* Create and initialize hashtable */
	dist_hashtable = create_hash_table(size * load_factor, &memory_heap, my_heap_size);
	upc_barrier;
	
	end_setup = UPC_TICKS_NOW();
   setupTime = UPC_TICKS_TO_SECS(end_setup - start_setup);

	if (MYTHREAD == 0)
	printf("Threads done with setting-up\n");
	
	/* Allocate local arrays used for book-keeping when using aggregated upc_memputs */
	list_t *local_buffs = (list_t*) malloc_chk(THREADS * CHUNK_SIZE * sizeof(list_t));

	UPC_INT64_T *local_index = (UPC_INT64_T*) malloc_chk(THREADS * sizeof(UPC_INT64_T));
	memset(local_index, 0, THREADS * sizeof(UPC_INT64_T));
	
	start_storing = UPC_TICKS_NOW();
		
	ptr = 0;
	idx = 0;
	list_t new_entry;

	UPC_INT64_T store_pos;
   char exts[2];
   
	while (ptr < kmers_read) {
		/* Set end of string chars to the appropriate positions:       */
		/* Pointer: working_buffer+ptr                ---> current kmer        */
		/* Pointer: working_buffer+ptr+KMER_LENGTH+1  ---> current extensions  */
		if (ufx_remote_thread[idx] >= 0) {
			// get cached remote thread
			remote_thread = ufx_remote_thread[idx];
			
#ifdef DEBUG
			/* if DEBUG, verify the calculated remote thread is actually the correct one */
#ifdef LSH
			hashval = LSHkmer(dist_hashtable->size, (char*) kmersarr[ptr]);
#else
			hashval = hashkmer(dist_hashtable->size, (char*) kmersarr[ptr]);
#endif
			assert( remote_thread == upc_threadof(&(dist_hashtable->table[hashval])) );
#endif
         packSequence((unsigned char*) (kmersarr[ptr]), new_entry.packed_key, KMER_LENGTH);
         
#ifdef MERACULOUS
         exts[0] = lefts[ptr];
         exts[1] = rights[ptr];
         new_entry.packed_extensions = convertExtensionsToPackedCode((unsigned char*) exts);
			new_entry.used_flag = UNUSED;
         new_entry.my_contig = NULL;
#endif
         
#ifdef MERDEPTH
         new_entry.count = counts[ptr];
#endif
         
#ifdef CEA
         new_entry.left_ext = lefts[ptr];
         new_entry.right_ext = rights[ptr];
#endif
			new_entry.next = NULL;
			
			/* Has more space at local buffer, so add to local buffer */
			if (local_index[remote_thread] <= CHUNK_SIZE - 1 ) {
				memcpy( &(local_buffs[local_index[remote_thread] + remote_thread * CHUNK_SIZE]), &new_entry, sizeof(list_t));
				local_index[remote_thread]++;
				if (VERBOSE > 2) LOG("Thread %d: sending %s to %d\n", MYTHREAD, kmersarr[ptr], remote_thread);
			}
			
			/* If buffer for that thread is full, do a remote upc_memput() */
			if (local_index[remote_thread] == CHUNK_SIZE) {
				store_pos = UPC_ATOMIC_FADD_I64((memory_heap.cached_heap_indices[remote_thread]), CHUNK_SIZE);
				upc_memput( (shared[] list_t*)  ((memory_heap.cached_heap_ptrs[remote_thread]) + store_pos)  , &(local_buffs[remote_thread * CHUNK_SIZE]), (CHUNK_SIZE) * sizeof(list_t));
				local_index[remote_thread] = 0;
			}
#ifdef DEBUG
			kmer_count++;
#endif
		} else {
			if (VERBOSE > 1) LOG("Thread %d: Skipping %.*s\n", MYTHREAD, KMER_LENGTH+3,kmersarr[ptr]);
		}

      ptr++;
		idx++;
	}
   
	assert(idx == my_ufx_lines);

	free(ufx_remote_thread);

	upc_barrier;
	
	/* Now check if there are few more kmers to be stored in the local buffs */
	for (i = 0; i < THREADS; i++) {
		/* Have to store local_index[i] items */
		if (local_index[i] != 0) {
			store_pos = UPC_ATOMIC_FADD_I64((memory_heap.cached_heap_indices[i]), local_index[i]);
			upc_memput( (shared[] list_t*)  ((memory_heap.cached_heap_ptrs[i]) + store_pos)  , &(local_buffs[i * CHUNK_SIZE]), (local_index[i]) * sizeof(list_t));
		}
	}
	
	upc_barrier;

	assert( memory_heap.heap_indices[MYTHREAD] == my_heap_size );
	/* Now each thread will iterate over its local heap and will add nodes to the local part of the hashtable - NO NEED FOR LOCKS!!! */
	int64_t heap_entry;
   
   list_t *local_filled_heap = (list_t*) memory_heap.cached_heap_ptrs[MYTHREAD];
   int64_t hashTableSize = dist_hashtable->size;
   shared[] list_t* shared_local_filled_heap = (memory_heap.heap_ptr[MYTHREAD]);
	//shared[] list_t* local_filled_heap = (memory_heap.heap_ptr[MYTHREAD]);
	unsigned char unpacked_kmer[KMER_LENGTH+1];
	unpacked_kmer[KMER_LENGTH] = '\0';

	for (heap_entry = 0; heap_entry < memory_heap.heap_indices[MYTHREAD]; heap_entry++) {
	
		unpackSequence((unsigned char*) &(local_filled_heap[heap_entry].packed_key), (unsigned char*) unpacked_kmer, KMER_LENGTH);
#ifdef LSH
      hashval = LSHkmer(dist_hashtable->size, (char*) (unpacked_kmer));
#else
      hashval = hashkmer(dist_hashtable->size, (char*) (unpacked_kmer));
#endif
		local_filled_heap[heap_entry].next = dist_hashtable->table[hashval].head;
		dist_hashtable->table[hashval].head =  (shared[] list_t*) &(shared_local_filled_heap[heap_entry]);

	}
	
	upc_barrier;

	end_storing = UPC_TICKS_NOW();
   storeTime = UPC_TICKS_TO_SECS(end_storing-start_storing);
	
#ifdef PROFILE	
	upc_barrier;
	
	if (MYTHREAD == 0) {
		printf("\n************* SET - UP TIME *****************");
		printf("\nTime spent on setting up the distributed hash table is %f seconds\n", UPC_TICKS_TO_SECS(end_setup -start_setup));

#ifdef DETAILED_IO_PROFILING
		printf("\n\n************* DETAILED TIMINGS FOR I/O AND STORING KMERS *****************\n");
#endif

	}
	
	upc_barrier;
	
#ifdef DETAILED_IO_PROFILING

	printf("Thread %d spent %f seconds on read I/O and %f seconds on storing %lld kmers\n", MYTHREAD, UPC_TICKS_TO_SECS(end_read-start_read), UPC_TICKS_TO_SECS(end_storing -start_storing));

#endif

	upc_barrier;
	
   DeAllocateAll(&kmersarr, &counts, &lefts, &rights, kmers_read);
   
#ifdef BUCKET_BALANCE_PROFILING
	if (MYTHREAD == 0) {
		printf("\n--------------------------------\n");
		printf("--------- Bucket balance -------\n");
		printf("--------------------------------\n");
	}
	if (MYTHREAD == 0) {
		for (i=0; i<THREADS; i++)
			printf("Thread %lld has %lld elements in its buckets\n", i, memory_heap.heap_indices[i]);		
	}
#endif

#endif

	upc_barrier;

	(*memory_heap_res) = memory_heap;
	return dist_hashtable;
}

#endif // BUILD_UFX_HASH_H

