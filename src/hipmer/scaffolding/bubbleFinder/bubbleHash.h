#ifndef BUBBLE_HASH_H
#define BUBBLE_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <upc.h>
#include <string.h>
#include <assert.h>

#include "bubbleFinder.h"


/* Creates and initializes a distributed hastable - collective function */
hash_table_t* create_hash_table(int64_t size, memory_heap_t *memory_heap, int64_t my_heap_size)
{
   hash_table_t *result;
   int64_t i;
   int64_t n_buckets = size;
   
   result = (hash_table_t*) malloc_chk(sizeof(hash_table_t));
   result->size = n_buckets;
   result->table = (shared[BS] bucket_t*) upc_all_alloc(n_buckets, sizeof(bucket_t));
   
   /* Initialization of every bucket */
   for (i = MYTHREAD; i < n_buckets; i+= THREADS) {
      result->table[i].head = NULL;
      result->table[i].bucket_lock = upc_global_lock_alloc();
   }
   
   memory_heap->my_heap_ptr = (shared[] list_t*) upc_alloc(my_heap_size * sizeof(list_t));
   memory_heap->my_data_ptr = (shared[] data_t*) upc_alloc(my_heap_size * sizeof(data_t));
   memory_heap->my_heap_pos = 0;
   memory_heap->my_data_pos = 0;
   
   upc_barrier;
   
   return result;
}

/* Creates pointsTo and initializes a distributed hastable - collective function */
phash_table_t* create_phash_table(int64_t size, pmemory_heap_t *memory_heap, int64_t my_heap_size)
{
   phash_table_t *result;
   int64_t i;
   int64_t n_buckets = size;
   
   result = (phash_table_t*) malloc_chk(sizeof(phash_table_t));
   result->size = n_buckets;
   result->table = (shared[BS] pbucket_t*) upc_all_alloc(n_buckets, sizeof(pbucket_t));
   
   /* Initialization of every bucket */
   for (i = MYTHREAD; i < n_buckets; i+= THREADS) {
      result->table[i].head = NULL;
      result->table[i].bucket_lock = upc_global_lock_alloc();
   }
   
   memory_heap->my_heap_ptr = (shared[] plist_t*) upc_alloc(my_heap_size * sizeof(plist_t));
   memory_heap->my_data_ptr = (shared[] pdata_t*) upc_alloc(my_heap_size * sizeof(pdata_t));
   memory_heap->my_heap_pos = 0;
   memory_heap->my_data_pos = 0;
   
   upc_barrier;
   
   return result;
}

/* Computes the hash value of null terminated sequence */
/*
int64_t hash(int64_t  hashtable_size, char *kmer)
{
    unsigned long hashval;
    hashval = 5381;
    for(; *kmer != '\0'; kmer++) hashval = (*kmer) +  (hashval << 5) + hashval;
      
    return hashval % hashtable_size;
}

int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
	unsigned long hashval;
	hashval = 5381;
	for(int i = 0; i < size; i++) {
		hashval = seq[i] +  (hashval << 5) + hashval;
	}
	
	return hashval % hashtable_size;
}
*/

/* Add entry to the hash table for tipInfo */
int add_tip(hash_table_t *hashtable, char *tip, int contig, memory_heap_t *memory_heap)
{
   char rcTip[2*KMER_LENGTH];
   char rcFlag = PLUS;
   char *tipKey = tip;
   reverseComplementSeq(tip, rcTip, 2*KMER_LENGTH);
   if ( memcmp(rcTip, tip, 2*KMER_LENGTH) < 0 ) {
      tipKey = rcTip;
      rcFlag = MINUS;
   }
   
   int64_t hashval = hashkey(hashtable->size, tipKey, 2*KMER_LENGTH);
   list_t cached_copy;
   bucket_t local_buc;
   shared[] list_t *result;
   
   /* Obtain the appropriate lock to ensure atomicity */
   upc_lock(hashtable->table[hashval].bucket_lock);
   
   local_buc = hashtable->table[hashval];
   result = local_buc.head;
   int found = 0;
   
   for (; result != NULL;) {
      cached_copy = *result;
      if (memcmp(cached_copy.tipKey, tipKey, 2*KMER_LENGTH) == 0) {
         found = 1;
         break;
      }
      result = cached_copy.next;
   }
   
   int64_t cur_data_pos = memory_heap->my_data_pos;
   shared[] data_t *my_data_ptr = memory_heap->my_data_ptr;
   int64_t cur_heap_pos;
   shared[] list_t *my_heap_ptr;
   
   if ( found == 1 ) {
      /* There is such a key, therefore just add the data in the list */
      my_data_ptr[cur_data_pos].next = cached_copy.data;
      my_data_ptr[cur_data_pos].contig = contig;
      my_data_ptr[cur_data_pos].rcFlag = rcFlag;
      result->nContigs++;
      result->data = &(my_data_ptr[cur_data_pos]);
      (memory_heap->my_data_pos)++;
   } else {
      /* There is not such a key, therefore add both the entry and data */
      cur_heap_pos = memory_heap->my_heap_pos;
      my_heap_ptr = memory_heap->my_heap_ptr;
      my_data_ptr[cur_data_pos].next = NULL;
      my_data_ptr[cur_data_pos].contig = contig;
      my_data_ptr[cur_data_pos].rcFlag = rcFlag;
      my_heap_ptr[cur_heap_pos].next = local_buc.head;
      memcpy( (char*) my_heap_ptr[cur_heap_pos].tipKey, tipKey, 2*KMER_LENGTH);
      my_heap_ptr[cur_heap_pos].nContigs = 1;
      my_heap_ptr[cur_heap_pos].data = &(my_data_ptr[cur_data_pos]);
      hashtable->table[hashval].head = &(my_heap_ptr[cur_heap_pos]);
      (memory_heap->my_heap_pos)++;
      (memory_heap->my_data_pos)++;
   }
   
   /* Now unlock the bucjet since it is safe to operate: The tip has been added successfully */
   upc_unlock(local_buc.bucket_lock);
   return 1;
}

shared[] list_t* lookup_and_copy_tipInfo(hash_table_t *hashtable, char *tipKey, list_t *copy)
{
   int64_t hashval = hashkey(hashtable->size, tipKey, 2*KMER_LENGTH);
   bucket_t local_buc;
   shared[] list_t *result;
   
   local_buc = hashtable->table[hashval];
   result = local_buc.head;
   
   for (; result != NULL;) {
      (*copy) = *result;
      if (memcmp(copy->tipKey, tipKey, 2*KMER_LENGTH) == 0) {
         return result;
      }
      result = copy->next;
   }
   
   return NULL;
}

shared[] plist_t* lookup_and_copy_pt(phash_table_t *hashtable, char *tipKey, plist_t *copy)
{
   int64_t hashval = hashkey(hashtable->size, tipKey, KMER_LENGTH);
   pbucket_t local_buc;
   shared[] plist_t *result;
   
   local_buc = hashtable->table[hashval];
   result = local_buc.head;
   
   for (; result != NULL;) {
      (*copy) = *result;
      if (memcmp(copy->tipKey, tipKey, KMER_LENGTH * sizeof(char)) == 0) {
         return result;
      }
      result = copy->next;
   }

   return NULL;
}

/* Add entry to the hash table for pointsTo */
int add_edge(phash_table_t *hashtable, char *tipKey, int objectID, char type, pmemory_heap_t *memory_heap)
{

   int64_t hashval = hashkey(hashtable->size, tipKey, KMER_LENGTH);
   plist_t cached_copy;
   pbucket_t local_buc;
   shared[] plist_t *result;
   
   /* Obtain the appropriate lock to ensure atomicity */
   upc_lock(hashtable->table[hashval].bucket_lock);
   
   local_buc = hashtable->table[hashval];
   result = local_buc.head;
   int found = 0;
   
   for (; result != NULL;) {
      cached_copy = *result;
      if (memcmp(cached_copy.tipKey, tipKey, KMER_LENGTH) == 0) {
         found = 1;
         break;
      }
      result = cached_copy.next;
   }
   
   int64_t cur_data_pos = memory_heap->my_data_pos;
   shared[] pdata_t *my_data_ptr = memory_heap->my_data_ptr;
   int64_t cur_heap_pos;
   shared[] plist_t *my_heap_ptr;
   
   if ( found == 1 ) {
      /* There is such a key, therefore just add the data in the list */
      my_data_ptr[cur_data_pos].next = cached_copy.data;
      my_data_ptr[cur_data_pos].object = objectID;
      my_data_ptr[cur_data_pos].type = type;
      result->nContigs++;
      result->data = &(my_data_ptr[cur_data_pos]);
      (memory_heap->my_data_pos)++;
   } else {
      /* There is not such a key, therefore add both the entry and data */
      cur_heap_pos = memory_heap->my_heap_pos;
      my_heap_ptr = memory_heap->my_heap_ptr;
      my_data_ptr[cur_data_pos].next = NULL;
      my_data_ptr[cur_data_pos].object = objectID;
      my_data_ptr[cur_data_pos].type = type;
      my_heap_ptr[cur_heap_pos].next = local_buc.head;
      memcpy( (char*) my_heap_ptr[cur_heap_pos].tipKey, tipKey, KMER_LENGTH * sizeof(char));
      my_heap_ptr[cur_heap_pos].nContigs = 1;
      my_heap_ptr[cur_heap_pos].data = &(my_data_ptr[cur_data_pos]);
      hashtable->table[hashval].head = &(my_heap_ptr[cur_heap_pos]);
      (memory_heap->my_heap_pos)++;
      (memory_heap->my_data_pos)++;
   }
   
   /* Now unlock the bucjet since it is safe to operate: The tip has been added successfully */
   upc_unlock(local_buc.bucket_lock);
   return 1;
}

int getNext(stream_t input, shared[1] extensions_t *linkertigs, shared[1] bubbletig_t *bubbletigs, phash_table_t *pointsTo, stream_t *output)
{
   int result;
   char type = input.type;
   char strand = input.strand;
   int64_t ID = input.objectId;
   char links[2*KMER_LENGTH];
   char in[KMER_LENGTH+1];
   char rcin[KMER_LENGTH+1];
   char out[KMER_LENGTH+1];
   char rcout[KMER_LENGTH+1];
   in[KMER_LENGTH] = '\0';
   rcin[KMER_LENGTH] = '\0';
   out[KMER_LENGTH] = '\0';
   rcout[KMER_LENGTH] = '\0';

   bubbletig_t cur_bubbletig;
   extensions_t cur_linkertig;
   plist_t local_pt;
   shared[] plist_t *res_ptr;
   pdata_t res_pdata;
   char pdata_type;
   
   result = 0;
   
   if (type == 'C') {
      cur_linkertig = linkertigs[ID];
      if (cur_linkertig.used_flag != USED_EXT) printf("FATAL ERROR - LINKERTIGS NOT DEFINED for [%ld]\n", ID);
      memcpy(in, &(cur_linkertig.data[0]), KMER_LENGTH * sizeof(char));
      memcpy(out, &(cur_linkertig.data[KMER_LENGTH]), KMER_LENGTH * sizeof(char));
   } else if (type == 'B') {
      cur_bubbletig = bubbletigs[ID];
      memcpy(in, &(cur_bubbletig.tipKey[0]), KMER_LENGTH * sizeof(char));
      memcpy(out, &(cur_bubbletig.tipKey[KMER_LENGTH]), KMER_LENGTH * sizeof(char));
   } else {
      printf("UNHANDLED CASE BUUUUUUUUUUUUUUUUUG\n");
   }
   
   reverseComplementSeq(in, rcout, KMER_LENGTH);
   reverseComplementSeq(out, rcin, KMER_LENGTH);
   res_ptr = NULL;
   
   if (strand == '-') {
      if ( (rcout[0] != '0') ) {
         res_ptr = lookup_and_copy_pt(pointsTo, rcout, &local_pt);
      } else {
         return 0;
      }
   } else {
      if ( (out[0] != '0') ) {
         res_ptr = lookup_and_copy_pt(pointsTo, out, &local_pt);
      } else {
         return 0;
      }
   }
   
   if (res_ptr == NULL) return 0;
   if ( local_pt.nContigs > 1 ) return 0;
   
   res_pdata = *(res_ptr->data);
   pdata_type = res_pdata.type;
   output->objectId = res_pdata.object;
   if (pdata_type == B_PLUS) {
      output->type = 'B';
      output->strand = '+';
   } else if (pdata_type == B_MINUS) {
      output->type = 'B';
      output->strand = '-';
   } else if (pdata_type == C_PLUS) {
      output->type = 'C';
      output->strand = '+';
   } else if (pdata_type == C_MINUS) {
      output->type = 'C';
      output->strand = '-';
   } else {
      printf("UNHANDLED CASE BUUUUUUUUUUUUUUUUUG\n");
   }
   
   return 1;
}


#endif // BUBBLE_HASH_H

