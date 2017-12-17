#ifndef __BUILD_SEED_INDEX_H
#define __BUILD_SEED_INDEX_H

#include "../common/common.h"
#include "../common/upc_compatibility.h"
extern int64_t seedsExtracted;
#define MAX_CONTIG_SIZE 900000

/* Builds the distributed seed index in parallel */
hash_table_t* buildSeedIndex(int64_t size, double exp_fact, char *inputContigsFileName, memory_heap_t *memory_heap_res, int subContigEffectiveLength, int64_t IDoffset ,int64_t myNSubcontigs, int read_length)
{
	hash_table_t *dist_hashtable;
	memory_heap_t memory_heap;
   FASTAFILE *ffp;
   char *contigBuffer;
   char *contigId;
   int contigID, contigLen;
   shared[] contig_t *cur_contig;
   char *contig;
   shared[] char *contig_content;
   int i, lex_ind, foundOccurence;
   char cur_kmer[KMER_LENGTH+1], auxiliary_kmer[KMER_LENGTH+1];
   char *kmer_to_store;
   auxiliary_kmer[KMER_LENGTH] = '\0';
   cur_kmer[KMER_LENGTH] = '\0';
   list_t new_entry;
	int64_t hashval;
	int64_t store_pos;
	int remote_thread, bufferOffset, startingSeedingPoint, endingSeedingPoint;
   contigInfo new_contig_info;
   shared[] list_t *curPtr;
   list_t curPtrCopy;
   int64_t subcontigsStoredSofar = 0;
   int64_t nCurSubContigs;
   int64_t x;
   int actualLength;
   char *contig_id;
	
	UPC_TICK_T start_read, end_read, start_storing, end_storing, start_setup, end_setup;
   start_setup = UPC_TICKS_NOW();
   
   /* Create and initialize hashtable */
	dist_hashtable = create_hash_table(size, &memory_heap, exp_fact);
   
	upc_barrier;
   
	end_setup = UPC_TICKS_NOW();

	if (MYTHREAD == 0)
		printf("Threads done with setting-up\n");
	
	/* Initialize lookup-table --> necessary for packing routines */
	upc_barrier;
	init_LookupTable();
	upc_barrier;
	  
	/* Allocate local arrays used for book-keeping when using aggregated upc_memputs */
	list_t *local_buffs = (list_t*) malloc_chk(THREADS * CHUNK_SIZE * sizeof(list_t));
	int64_t *local_index = (int64_t*) malloc_chk(THREADS * sizeof(int64_t));
	memset(local_index, 0, THREADS * sizeof(int64_t));
	
	start_storing = UPC_TICKS_NOW();

   ffp = OpenFASTA(inputContigsFileName);
   subcontigsStoredSofar = 0;
   double overlap = 2 * (read_length + SLACK) - KMER_LENGTH;

   if (!subContigEffectiveLength)
       printf("Skipping chopping optimization\n");

   while ( subcontigsStoredSofar != myNSubcontigs && ReadFASTA(ffp, &contigBuffer, &contig_id, &contigLen) ) {
      /* Process a Contig */
      contigID = atoi(contig_id+7);      // Assume target names are in the form ContigXXXXXX
      
      /* Chop the contig into smaller subpieces */
      if (!subContigEffectiveLength)
          nCurSubContigs = 1;
      else
          nCurSubContigs = (int64_t) ceil((double) contigLen/ (double) subContigEffectiveLength);
      
      for ( x = 0; x < nCurSubContigs; x++) {
          if (nCurSubContigs == 1) {
              actualLength = contigLen;
          } else {
              /* The first  subcontig of a contig does not suffer from the "overlap" overhead */
              /* Note: OVERLAP should be 2*(read_length+slack) - KMER_LENGTH */
              actualLength = (x == 0) ? subContigEffectiveLength : subContigEffectiveLength + overlap;
              if (x == (nCurSubContigs - 1)) actualLength = contigLen - subContigEffectiveLength * (nCurSubContigs-1) + overlap;
          }
         
         /* Allocate space for the contig in shared memory */
         cur_contig = (shared[] contig_t*) upc_alloc(sizeof(contig_t));
         contig_content = (shared[] char*) upc_alloc((actualLength+1) * sizeof(char));
         bufferOffset = (x == 0) ? 0 : x * subContigEffectiveLength - overlap;
         memcpy((char*)contig_content, contigBuffer + bufferOffset, actualLength * sizeof(char));
         contig_content[actualLength] = '\0';
         cur_contig->contig = contig_content;
         cur_contig->contig_id = subcontigsStoredSofar + IDoffset;
         subcontigsStoredSofar++;
         cur_contig->length = actualLength;
         cur_contig->parentLength = contigLen;
         cur_contig->parentContigID = contigID;
         cur_contig->offsetInParent = bufferOffset;
         cur_contig->uuType = UU_CONTIG;
         contig = (char*) contig_content;
         
         startingSeedingPoint = (x == 0) ? 0 : read_length - KMER_LENGTH + SLACK;
         endingSeedingPoint = actualLength - 1 - (read_length + SLACK);
         if (x == (nCurSubContigs-1)) endingSeedingPoint = actualLength - KMER_LENGTH;

         /* Mark the k-mers of the contig appropriately */
         for (i=startingSeedingPoint; i <= endingSeedingPoint; i++) {
            memcpy(cur_kmer, contig+i, KMER_LENGTH * sizeof(char));
            
            /* Search for the canonical kmer */
            reverseComplementKmer(cur_kmer, auxiliary_kmer);
            lex_ind = strcmp(cur_kmer, auxiliary_kmer);
            if ( lex_ind < 0)
               kmer_to_store = cur_kmer;
            else  {
               kmer_to_store = auxiliary_kmer;
            }
            
            /* Pack and insert kmer info into local_list_t buffer */
            hashval = hashstr(dist_hashtable->size, (char*)kmer_to_store);
            remote_thread = upc_threadof(&(dist_hashtable->table[hashval]));
            packSequence((unsigned char*) kmer_to_store, new_entry.packed_key, KMER_LENGTH);
            new_entry.next = NULL;
            new_entry.contigInfoArrayPtr = NULL;
            new_contig_info.posInContig = i;
            new_contig_info.my_contig = cur_contig;
            new_entry.firstContig = new_contig_info;
            new_entry.nExtraContigs = 0;
            new_entry.posInSubArray = 0;
            
            seedsExtracted++;
            
            /* Has more space at local buffer, so add to local buffer */
            if (local_index[remote_thread] <= CHUNK_SIZE - 1 ) {
               memcpy( &(local_buffs[local_index[remote_thread] + remote_thread * CHUNK_SIZE]), &new_entry, sizeof(list_t));
               local_index[remote_thread]++;
            }
            
            /* If buffer for that thread is full, do a remote upc_memput() */
            if (local_index[remote_thread] == CHUNK_SIZE) {
               store_pos = bupc_atomicI64_fetchadd_strict(&(memory_heap.heap_indices[remote_thread]), CHUNK_SIZE);
               assert(store_pos + CHUNK_SIZE <= memory_heap.heap_ptr_sizes[remote_thread]);
               upc_memput( (shared[] list_t*)  ((memory_heap.heap_ptr[remote_thread]) + store_pos)  , &(local_buffs[remote_thread * CHUNK_SIZE]), (CHUNK_SIZE) * sizeof(list_t));
               local_index[remote_thread] = 0;
            }
         }
      }
      
      free(contigBuffer);
      free(contig_id);
   }
   
   if(ffp!=NULL)
       CloseFASTA(ffp);
   //printf("Thread %d done with seeding part\n", MYTHREAD);
   /* Sanity check */
   if (subcontigsStoredSofar != myNSubcontigs) printf("Fatal ERROR: Incorrectly choped parent contig %ld VS %ld !\n", subcontigsStoredSofar, myNSubcontigs);
   
	upc_barrier;
   if (MYTHREAD == 0) printf("Done allocating and reading contigs\n");
	
	/* Now check if there are a few more kmers to be stored in the local buffs */
	for (i = 0; i < THREADS; i++) {
		/* Have to store local_index[i] items */
		if (local_index[i] != 0) {
			store_pos = bupc_atomicI64_fetchadd_strict(&(memory_heap.heap_indices[i]), local_index[i]);
                        assert(store_pos + local_index[i] <= memory_heap.heap_ptr_sizes[i]);
			upc_memput( (shared[] list_t*)  ((memory_heap.heap_ptr[i]) + store_pos)  , &(local_buffs[i * CHUNK_SIZE]), (local_index[i]) * sizeof(list_t));
		}
	}
	
	upc_barrier;

        if (memory_heap.heap_indices[MYTHREAD] > memory_heap.heap_ptr_sizes[MYTHREAD]) {
            fprintf(stderr, "Memory over run on thread %d. Increase the expansion factor for the meraligner dist hashtable!\n", MYTHREAD);
            upc_global_exit(1);
        }
   
   /* First pass over the seed heap to insert then first occurence of each seed in the hash-table and to count the number of occurences */

   int64_t heap_entry;
   int64_t recvMultipleLocs = 0, sendMultipleLocs = 0;
	
	shared[] list_t* local_filled_heap = (memory_heap.heap_ptr[MYTHREAD]);
        assert( upc_threadof(local_filled_heap) == MYTHREAD );
	unsigned char unpacked_kmer[KMER_LENGTH+1];
	unpacked_kmer[KMER_LENGTH] = '\0';
   shared[] contigInfo *extraContigInfoArray = NULL;
   
   for (heap_entry = 0; heap_entry < memory_heap.heap_indices[MYTHREAD]; heap_entry++) {
	unpackSequence((unsigned char*) &(local_filled_heap[heap_entry].packed_key), (unsigned char*) unpacked_kmer, KMER_LENGTH);
	hashval = hashstr(dist_hashtable->size, (char*)(unpacked_kmer));
      /* Check if already an occurence of the seed exists in the hash table */
      assert( upc_threadof( &(dist_hashtable->table[hashval]) ) == MYTHREAD );
      curPtr = dist_hashtable->table[hashval].head;
      foundOccurence = FALSE;
      while ( (foundOccurence == FALSE) && (curPtr !=NULL) ) {
         curPtrCopy = *curPtr;
         if (memcmp( (unsigned char*) local_filled_heap[heap_entry].packed_key, (unsigned char*) curPtrCopy.packed_key, KMER_PACKED_LENGTH * sizeof(unsigned char) ) == 0 ) {
            assert( curPtr != &(local_filled_heap[heap_entry]) );
            assert( curPtr->nExtraContigs >= 0 );
            assert( local_filled_heap[heap_entry].nExtraContigs == 0 );

            foundOccurence = TRUE;
#ifdef STORE_OPT
            assert( upc_threadof(curPtr) == MYTHREAD );
            int oldExtra = (curPtr->nExtraContigs)++;
#else
            int oldExtra = UPC_ATOMIC_FADD_I64(&(curPtr->nExtraContigs), 1);
#endif
            assert(oldExtra + 1 == curPtr->nExtraContigs);
            assert(local_filled_heap[heap_entry].nExtraContigs == 0);
            local_filled_heap[heap_entry].nExtraContigs = -1; // Mark that this entry contains additional info!!!
            sendMultipleLocs++;
         } else {
            curPtr = curPtrCopy.next;
         }
      }
      
      /* If not already in the hashtable, add it there */
      if (foundOccurence == FALSE) {
         assert( upc_threadof(&(dist_hashtable->table[hashval])) == MYTHREAD );
         assert( upc_threadof(&local_filled_heap[heap_entry]) == MYTHREAD );
         local_filled_heap[heap_entry].next = dist_hashtable->table[hashval].head;
         dist_hashtable->table[hashval].head = &local_filled_heap[heap_entry];
      }
   }

   /* Second pass over the seed heap to count the number of extraContigs to allocate */
   for (heap_entry = 0; heap_entry < memory_heap.heap_indices[MYTHREAD]; heap_entry++) {
      assert( upc_threadof( &local_filled_heap[heap_entry] ) == MYTHREAD );
      if (local_filled_heap[heap_entry].nExtraContigs > 0) {
          recvMultipleLocs += local_filled_heap[heap_entry].nExtraContigs;
      }
   }
   assert(recvMultipleLocs == sendMultipleLocs); /* all operations are local */
#ifdef DEBUG
   fprintf(stderr, "Thread %d: sending %lld, receiving %lld\n", MYTHREAD, (long long) sendMultipleLocs, (long long) recvMultipleLocs); 
#endif
   
   /* Allocate an array for the extra contigs */
   assert(memory_heap.contigInfoArrays[MYTHREAD] == NULL);
   if (recvMultipleLocs != 0) {
      extraContigInfoArray = (shared[] contigInfo*) upc_alloc(recvMultipleLocs * sizeof(contigInfo));
      memory_heap.contigInfoArrays[MYTHREAD] = extraContigInfoArray;
   }
   int64_t posInBuffer = 0;
   //printf("Thread ID is %d identified %lld multiple locations\n", MYTHREAD, nMultipleLocs);

   /* Third pass over the seed heap to allocate appropriately the information for the contigInfoArray from allocation */
   for (heap_entry = 0; heap_entry < memory_heap.heap_indices[MYTHREAD]; heap_entry++) {
      /* Should do some action only if extra contigs exist and this is the instance in the table */
      if (local_filled_heap[heap_entry].nExtraContigs > 0) {
         assert(upc_threadof(local_filled_heap + heap_entry) == MYTHREAD);
         assert(posInBuffer + local_filled_heap[heap_entry].nExtraContigs <= recvMultipleLocs);
         /* We should do the allocation - This is the first occurence of the seed */
         assert(extraContigInfoArray != NULL);
         assert( upc_threadof( &(extraContigInfoArray[posInBuffer]) ) == MYTHREAD );
         local_filled_heap[heap_entry].contigInfoArrayPtr = (shared[] contigInfo*) &(extraContigInfoArray[posInBuffer]);
         posInBuffer += local_filled_heap[heap_entry].nExtraContigs;
         assert(local_filled_heap[heap_entry].posInSubArray == 0);
      }
   }
   assert(posInBuffer == recvMultipleLocs);
   int64_t sent = 0;

   /* Fourth pass over the seed heap to insert appropriately the information for the contigInfoArray */
   for (heap_entry = 0; heap_entry < memory_heap.heap_indices[MYTHREAD]; heap_entry++) {
      
      if (local_filled_heap[heap_entry].nExtraContigs == -1) {
         /* This includes extra contig info, so find the appropriate entry in the hash table and add that info */
         unpackSequence((unsigned char*) &(local_filled_heap[heap_entry].packed_key), (unsigned char*) unpacked_kmer, KMER_LENGTH);
         hashval = hashstr(dist_hashtable->size, (char*)(unpacked_kmer));
         curPtr = dist_hashtable->table[hashval].head;
         foundOccurence = FALSE;
         while ( (foundOccurence == FALSE) && (curPtr !=NULL) ) {
            curPtrCopy = *curPtr;
            assert( upc_threadof( local_filled_heap[heap_entry].packed_key ) == MYTHREAD );
            if (memcmp( (unsigned char*) local_filled_heap[heap_entry].packed_key, (unsigned char*) curPtrCopy.packed_key, KMER_PACKED_LENGTH * sizeof(unsigned char) ) == 0 ) {
               assert(curPtrCopy.contigInfoArrayPtr != NULL);
               assert(curPtrCopy.contigInfoArrayPtr == curPtr->contigInfoArrayPtr);
               foundOccurence = TRUE;
#ifdef STORE_OPT
               assert(upc_threadof(curPtr) == MYTHREAD);
               int posInSubArray = (curPtr->posInSubArray)++;
#else
               int posInSubArray = UPC_ATOMIC_FADD_I64(&(curPtr->posInSubArray), 1);
#endif
               assert(posInSubArray == curPtrCopy.posInSubArray);
               assert(posInSubArray + 1 == curPtr->posInSubArray);
               assert(posInSubArray < curPtrCopy.nExtraContigs);
               assert(upc_threadof( &(curPtrCopy.contigInfoArrayPtr[posInSubArray]) ) == MYTHREAD);
               assert(upc_threadof( &(local_filled_heap[heap_entry].firstContig) ) == MYTHREAD);
               //curPtrCopy.contigInfoArrayPtr[posInSubArray] = local_filled_heap[heap_entry].firstContig;
               memcpy( (void *) &(curPtrCopy.contigInfoArrayPtr[posInSubArray]), (void*) &local_filled_heap[heap_entry].firstContig, sizeof(contigInfo) );
               sent++;
            } else {
               curPtr = curPtrCopy.next;
            }
         }
         assert(foundOccurence);
      }
      /* If nExtraContigs != 0 then mark the corresponding contig as NON_UU_CONTIG */
      if (local_filled_heap[heap_entry].nExtraContigs != 0) {
         bupc_atomicI64_cswap_strict(&(local_filled_heap[heap_entry].firstContig.my_contig->uuType), UU_CONTIG, NON_UU_CONTIG);
      }
   }
   assert(sent == sendMultipleLocs);
   upc_fence;
   end_storing = UPC_TICKS_NOW();
	
#ifdef PROFILE
	upc_barrier;
	
	if (MYTHREAD == 0) {
		printf("\n************* SET - UP TIME *****************");
		printf("\nTime spent on setting up the distributed hash table is %d seconds\n", ((int)UPC_TICKS_TO_SECS(end_setup -start_setup)));
	}
#endif
   
	upc_barrier;
   
	(*memory_heap_res) = memory_heap;
	return dist_hashtable;
}

#endif
