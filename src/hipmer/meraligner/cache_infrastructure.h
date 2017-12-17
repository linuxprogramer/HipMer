#ifndef __CACHE_INFRASTRUCTURE_H
#define  __CACHE_INFRASTRUCTURE_H

#include "../common/common.h"

#define UNUSED_SLOT 0
#define USED_SLOT 1
#define UNREADY_SLOT 2

typedef shared[] char* contigDataPtr;
typedef shared[] list_t kmerData;
typedef shared[] list_t* kmerDataPtr;

/* Computes the hash value of a k-mer key for the cache */
/*
int64_t cacheHash(int64_t cacheCapacity, char *kmer)
{
   unsigned long hashval;
   hashval = 5381;
   for(; *kmer != '\0'; kmer++) hashval = (*kmer) +  (hashval << 5) + hashval;
   
   return hashval % cacheCapacity;
}
*/


/* Collective function that creates the cache infrastructure */
int createCache(int64_t nContigsIn, shared int64_t **cacheSize, shared[1] contigDataPtr **cacheTable, int64_t capacityPerCache)
{
   int64_t nContigs = (nContigsIn / cores_per_node) * cores_per_node + cores_per_node ;

   /* Check if caches have been already created */
   if ((*cacheTable) != NULL) {
      return 1;
   }
   
   shared[1] int64_t *cacheSizesArray;
   shared[1] contigDataPtr *tmpCacheTable;
   int64_t i;
   int64_t myNode = MYTHREAD / cores_per_node;
   int64_t nodes = (int64_t) ceil((double) THREADS / (double)cores_per_node);
   int64_t cacheTablesSize = nodes * nContigs;
   
   /* Create cache tables that will contain pointers to node-local contigs */
   if (MYTHREAD == 0)
      printf("Cache size meta-data is %f Gigabytes\n", (cacheTablesSize * sizeof(contigDataPtr))/1024.0/1024.0/1024.0);
   tmpCacheTable = (shared[1] contigDataPtr*) upc_all_alloc(cacheTablesSize , sizeof(contigDataPtr));
   if (tmpCacheTable == NULL) 
       DIE("Failed to allocate %ld of %lu for cacheTable\n", cacheTablesSize, (unsigned long)sizeof(contigDataPtr));

   for (i=MYTHREAD; i < cacheTablesSize; i+=THREADS)
      tmpCacheTable[i] = NULL;
   (*cacheTable) = tmpCacheTable;
   
   /* Create array with int values that indicate the current size of node-local cache */
   cacheSizesArray = (shared[1] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
   if (tmpCacheTable == NULL) 
       DIE("Failed to allocate cache sizes\n");
   cacheSizesArray[MYTHREAD] = capacityPerCache;
   upc_barrier;
   (*cacheSize) = &cacheSizesArray[myNode * cores_per_node];
   cache_hits = 0;
   cache_queries = 0;
   
   return 1;
}

contigDataPtr findRemoteContigPtr(int contigRequested, int contigLength, shared[] contig_t *contig_ptr, int start, shared int64_t *cachePtr, shared[1] contigDataPtr *cacheTable)
{
   contigDataPtr resultPtr;
   
   /* Look up node-local cache */
   int64_t chunk = contigRequested / cores_per_node;
   int64_t offset = contigRequested % cores_per_node;
   int64_t nid = MYTHREAD / cores_per_node;
   int64_t nodes = (int64_t) ceil((double) THREADS / (double)cores_per_node);
   int64_t loc = nid * cores_per_node + chunk * THREADS + offset;
   int64_t spaceInCache;
   shared[] char *tmp_res;
   
   //cache_queries++;
   resultPtr = cacheTable[loc];
   if (resultPtr == NULL) {
      /* Cache miss: Atomically check if there is sufficient size in the node-local cache */
      spaceInCache = bupc_atomicI64_fetchadd_strict(cachePtr, (-1) * contigLength);
      if ( spaceInCache >= contigLength ) {
         tmp_res = (shared[] char*) upc_alloc(contigLength * sizeof(char));
         upc_memget((char*)tmp_res, &(contig_ptr->contig[start]), contigLength * sizeof(char));
         cacheTable[loc] = (contigDataPtr) tmp_res;
         resultPtr = (contigDataPtr) tmp_res;
         upc_fence;
      }
      else {
         /* Since we did not use the space, add it back  */
         spaceInCache = bupc_atomicI64_fetchadd_strict(cachePtr, contigLength);
         resultPtr = (contigDataPtr) (&(contig_ptr->contig[start]));
      }
   } else {
      /* Cache hit */
      //cache_hits++;
   }
   
   return resultPtr;
}

/* Collective function that creates the cache infrastructure for k-mers */
/* Returns the number of slots for k-mers per node cache */
int64_t createKmerCache(shared[1] list_t **cacheTable, int64_t capacityPerCache, shared[1] int64_t **cacheFlags)
{
   int64_t datasize = sizeof(list_t) + sizeof(int64_t);
   int64_t nSlots = (int64_t) floor((double) capacityPerCache / (double) datasize );
   nSlots = (nSlots / cores_per_node) * cores_per_node + cores_per_node;
   int64_t nSlotsPerthread = (int64_t) floor((double) nSlots/ (double) cores_per_node );
   shared[1] list_t *tmpCacheTable;
   shared[1] int64_t *tmpCacheFlags;

   int64_t i;
   int64_t myNode = MYTHREAD / cores_per_node;
   int64_t nodes = (int64_t) ceil((double) THREADS / (double)cores_per_node);
   int64_t cacheTablesSize = nodes * nSlots;
   
   /* Create cache tables that will contain pointers to node-local contigs */
   if (MYTHREAD == 0)
      printf("Cache table and cache flags is %f Gigabytes (%f / node)\n", (cacheTablesSize * datasize)/1024.0/1024.0/1024.0, (cacheTablesSize * datasize)/1024.0/1024.0/1024.0/(THREADS/cores_per_node));
      
   tmpCacheTable = (shared[1] list_t*) upc_all_alloc(cacheTablesSize , sizeof(list_t));
   tmpCacheFlags = (shared[1] int64_t*) upc_all_alloc(cacheTablesSize , sizeof(int64_t));
   if (tmpCacheTable == NULL || tmpCacheFlags == NULL) 
       DIE("Failed to allocate %ld kmer cache table entries\n", cacheTablesSize); 
   
   /* Initialize with guard values */
   for (i = 0; i < nSlotsPerthread; i++) {
      tmpCacheFlags[MYTHREAD + i * THREADS] = UNUSED_SLOT;
   }
   
   upc_barrier;

   if (MYTHREAD == 0)
       printf("Cache initialized\n");
   
   (*cacheTable) = tmpCacheTable;
   (*cacheFlags) = tmpCacheFlags;
   
   return nSlots;
}

/* Returns a pointer to the requested k-mer. If the k-mer is cached then the pinter is node-local. If not, then the processor tries to cache the kmer in the node-local cache. Eventually the pointer is a valid one shared ptr to the requested k-mer (local or not)  */
kmerDataPtr lookupKmerInCache(hash_table_t *hashtable, unsigned char *kmer, int64_t cacheCapacity, shared[1] list_t *cacheTable, shared[1] int64_t *cacheFlags)
{
   
   /* Look up node-local cache */
   int64_t kmerRequested = hashstr(cacheCapacity, (void*)kmer);
   int64_t chunk = kmerRequested / cores_per_node;
   int64_t offset = kmerRequested % cores_per_node;
   int64_t nid = MYTHREAD / cores_per_node;
   int64_t nodes = (int64_t) ceil((double) THREADS / (double)cores_per_node);
   int64_t loc = nid * cores_per_node + chunk * THREADS + offset;
   unsigned char packed_key[KMER_PACKED_LENGTH];
   kmerDataPtr remoteKmerPtr;
   kmerDataPtr resultPtr;
   list_t localKmer;
   int64_t slotState;
   
   cache_queries++;
   
   slotState = bupc_atomicI64_cswap_strict(&cacheFlags[loc], UNUSED_SLOT, UNREADY_SLOT);
   if (slotState == UNUSED_SLOT) {
      /* In this case the k-mer is not in the cache and the cache-slot is empty -- look it up in the global hash table */
      remoteKmerPtr = lookup_kmer(hashtable, kmer);
      if (remoteKmerPtr == NULL) {
         /* The k-mer was not found in the global hash table, so return NULL     */
         /* Mark the slot ---ATOMICALLY--- as unused since no k-mer was found... */
         slotState = bupc_atomicI64_cswap_strict(&cacheFlags[loc], UNREADY_SLOT, UNUSED_SLOT);
         resultPtr = NULL;
      } else {
         /* Store the hash entry in the cache-slot */
         cacheTable[loc] = *remoteKmerPtr;
         resultPtr = (kmerDataPtr) &cacheTable[loc];
         slotState = bupc_atomicI64_cswap_strict(&cacheFlags[loc], UNREADY_SLOT, USED_SLOT);
      }
   } else if (slotState == UNREADY_SLOT) {
      /* Another processor just claimed that slot, so do a remote look-up ... */
      resultPtr = lookup_kmer(hashtable, kmer);
   } else if (slotState == USED_SLOT) {
      /* The cache-slot is not empty - check if this is the requested kmer */
      localKmer = cacheTable[loc];
      packSequence(kmer, packed_key, KMER_LENGTH);
      if (comparePackedSeq(packed_key, localKmer.packed_key, KMER_PACKED_LENGTH) == 0) {
         /* This is indeed the requested k-mer */
			resultPtr = (kmerDataPtr) &cacheTable[loc];
         cache_hits++;
		} else {
         /* The slot is occupied by another k-mer than the one requested... */
         resultPtr = lookup_kmer(hashtable, kmer);
      }
   }
   
   return resultPtr;
}

#endif
