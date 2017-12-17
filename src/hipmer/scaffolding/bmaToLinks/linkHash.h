#define SPLINT 0
#define SPAN 1
#define LINK_CHUNK_SIZE 100
#define BS 1
#define PLUS 0
#define MINUS 1
#include <upc.h>
#include "bma_meta.h"
#include "../../common/common.h"

typedef struct link_t link_t;
struct link_t{
   int end1_id;                  // ID of link's endpoint 1
   int end2_id;                  // ID of link's endpoint 2
   char compressedTails;         // Code that repesents compressed info for (*.3 / *.5)
   char nature;
   unsigned char link_type;      // Type of link: SPLINT or SPAN
   int lib_id;                   // Library id
   int endSeparation;
   int d1;
   int d2;
   int o1_length;
   int o2_length;
};

typedef shared[] link_t* shared_link_ptr;

/* Link heap data structure */
typedef struct link_heap_t link_heap_t;
struct link_heap_t {
	shared[BS] shared_link_ptr *link_ptr;					// Pointers to shared memory heaps
	shared[BS] int64_t *heap_indices;						// Indices of remote heaps
};

/* Creates shared heaps that required for links hash table - IMPORTANT: This is a collective function */
int create_link_heaps(int64_t size, link_heap_t *link_heap)
{
#ifdef DEBUG
        printf("Thread %d: Allocating %ld link heaps\n", MYTHREAD, size);
#endif
	link_heap->link_ptr = (shared[BS] shared_link_ptr*) upc_all_alloc(THREADS, sizeof(shared_link_ptr));
    link_heap->heap_indices = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
	link_heap->link_ptr[MYTHREAD] = (shared[] link_t*) upc_alloc( (size + THREADS) / THREADS * sizeof(link_t));
    if (link_heap->link_ptr == NULL || link_heap->link_ptr[MYTHREAD] == NULL || link_heap->heap_indices == NULL) 
        DIE("Could not allocate link_heaps!\n");
	link_heap->heap_indices[MYTHREAD] = 0;
	upc_barrier;
	
	return 0;
}

/* Allocate local arrays used for book-keeping when using aggregated upc_memputs */
int allocate_local_buffs(link_t **local_buffs, int64_t **local_index)
{
	(*local_buffs) = (link_t*) malloc_chk(THREADS * LINK_CHUNK_SIZE * sizeof(link_t));
	(*local_index) = (int64_t*) malloc_chk(THREADS * sizeof(int64_t));
	memset((*local_index), 0, THREADS * sizeof(int64_t));
   return 0;
}

/* Adds a link to the shared link heap */
int add_link_to_shared_heaps(link_t *new_entry, int64_t hashval, int64_t *local_index, link_t *local_buffs, link_heap_t link_heap)
{
   int64_t store_pos;
   int remote_thread = hashval % THREADS;
   
   /* Store link first to local buffer designated for remote thread */
   if (local_index[remote_thread] <= LINK_CHUNK_SIZE - 1) {
      memcpy( &(local_buffs[local_index[remote_thread] + remote_thread * LINK_CHUNK_SIZE]), new_entry, sizeof(link_t));
      local_index[remote_thread]++;
   }
   
   /* If buffer for that thread is full, do a remote upc_memput() */
   if (local_index[remote_thread] == LINK_CHUNK_SIZE) {
      store_pos = bupc_atomicI64_fetchadd_strict(&(link_heap.heap_indices[remote_thread]), LINK_CHUNK_SIZE);
      upc_memput( (shared[] link_t*)  ((link_heap.link_ptr[remote_thread]) + store_pos)  , &(local_buffs[remote_thread * LINK_CHUNK_SIZE]), (LINK_CHUNK_SIZE) * sizeof(link_t));
      local_index[remote_thread] = 0;
   }

   return 0;
}

/* Adds remaining links to the shared heaps. Should be called when all calls "add_link_to_shared_heaps()" have been done */
int add_rest_links_to_shared_heaps(int64_t *local_index, link_t *local_buffs, link_heap_t link_heap)
{
   int i;
   int64_t store_pos;
   
   for (i = 0; i < THREADS; i++) {
		if (local_index[i] != 0) {
			store_pos = bupc_atomicI64_fetchadd_strict(&(link_heap.heap_indices[i]), local_index[i]);
			upc_memput( (shared[] link_t*)  ((link_heap.link_ptr[i]) + store_pos)  , &(local_buffs[i * LINK_CHUNK_SIZE]), (local_index[i]) * sizeof(link_t));
		}
	}
   
   return 0;
}
/*
int64_t linkhash(int64_t linkhash_size, char *linkname)
{
   unsigned long hashval;
   hashval = 5381;
   for(; *linkname != '\0'; linkname++) hashval = (*linkname) +  (hashval << 5) + hashval;
   return hashval % linkhash_size;
}
*/

/* Data structure to store data regarding a link */
typedef struct datum_t datum_t;
struct datum_t{
   int endSeparation;            // End Separation field of dataum
   int d1;                       // Span specific data
   int d2;                       // Span specific data
   int lib_id;                   // Library identifier
   int count;
	datum_t *next;                // Pointer to next entry in the same datum list
};

/* Essentially the key now is the combination: end1_id | end2_id | linktype */
typedef struct linkList_t linkList_t;
struct linkList_t{
   int end1_id;                  // ID of link's endpoint 1
   int end2_id;                  // ID of link's endpoint 2
   char compressedTails;
   unsigned char link_type;      // Type of link: SPLINT or SPAN
   int o1_length;
   int o2_length;
   datum_t *data;                // Data list
	linkList_t *next;					// Pointer to next entry in the same bucket
};

typedef struct link_hash_table_t_ {
   int64_t size;                     /* the size of the table */
   linkList_t **table;           /* the table elements */
} link_hash_table_t;

/* Create a link_hash_table */
link_hash_table_t *create_link_hash_table(int64_t size)
{
   link_hash_table_t *new_table;
   int64_t i;
   
   if (size<1) size = 1;
   
   /* Attempt to allocate memory for the table structure */
   if ((new_table = (link_hash_table_t*) malloc_chk(sizeof(link_hash_table_t))) == NULL) {
      return NULL;
   }
   
   /* Attempt to allocate memory for the table itself */
   if ((new_table->table = malloc_chk(sizeof(linkList_t*) * size)) == NULL) {
      return NULL;
   }
   
   /* Initialize the elements of the table */
   for(i=0; i<size; i++) new_table->table[i] = NULL;
   
   /* Set the table's size */
   new_table->size = size;
   
   return new_table;
}

char compressTails(char tail1, char tail2) {
   char result = 0;
   
   if ((tail1 == '3') && (tail2 == '3')) {
      result = 0;
   }
   
   if ((tail1 == '3') && (tail2 == '5')) {
      result = 1;
   }
   
   if ((tail1 == '5') && (tail2 == '3')) {
      result = 2;
   }
   
   if ((tail1 == '5') && (tail2 == '5')) {
      result = 3;
   }
   
   return result;
}

int64_t link_hash_val(int e1, int e2, char compressedTails, int64_t linkhash_size)
{
   char linkname[LINK_CHUNK_SIZE];
   int tail1, tail2;
   if ( compressedTails == 0 ) {
      tail1 =3;
      tail2 =3;
   }
   if ( compressedTails == 1 ) {
      tail1 =3;
      tail2 =5;
   }
   if ( compressedTails == 2 ) {
      tail1 =5;
      tail2 =3;
   }
   if ( compressedTails == 3 ) {
      tail1 =5;
      tail2 =5;
   }
   
   char linkAux1[LINK_CHUNK_SIZE];
   char linkAux2[LINK_CHUNK_SIZE];
   
   sprintf_chk(linkAux1, LINK_CHUNK_SIZE, "Contig%d.%d", e1, tail1);
   sprintf_chk(linkAux2, LINK_CHUNK_SIZE, "Contig%d.%d", e2, tail2);
   
   if (strcmp(linkAux1,  linkAux2) < 0) {
      sprintf_chk(linkname, LINK_CHUNK_SIZE, "Contig%d.%d<=>Contig%d.%d", e1, tail1, e2, tail2);
   }
   else {
      sprintf_chk(linkname, LINK_CHUNK_SIZE, "Contig%d.%d<=>Contig%d.%d", e2, tail2, e1, tail1);
   }
   //printf("Example linkname %s\n", linkname);
   
   
   return hashstr(linkhash_size, linkname);
}

/* Lookup function for link hash table */
linkList_t *lookup_link(link_hash_table_t *hashtable, int e1, int e2, char compressedTails ,unsigned char type, int64_t hashval)
{
   linkList_t *list;
   
   for(list = hashtable->table[hashval]; list != NULL; list = list->next) {
      if ((e1 == list->end1_id) && (e2 == list->end2_id) && (type ==list->link_type) && (compressedTails == list->compressedTails)) return list;
   }
   return NULL;
}

/* Implement efficient memory allocator to avoid multiple malloc() calls */
int create_local_heaps(int64_t max_size, datum_t **datum_heap, linkList_t **linkList_heap)
{
   if(max_size<1) max_size = 1;
   (*datum_heap) = (datum_t*) malloc_chk(max_size * sizeof(datum_t));
   (*linkList_heap) = (linkList_t*) malloc_chk(max_size * sizeof(linkList_t));
   return 0;
}

/* Add link data to the hashtable */
int add_link_data(link_hash_table_t *hashtable, link_t *cur_link, int64_t *datum_heap_pos, datum_t *datum_heap, int64_t *link_heap_pos, linkList_t *linkList_heap, bma_info *bmas, int64_t *nLinks, FILE *myLog)
{
   linkList_t *new_bucket;
   datum_t *new_datum, *span_data;
   linkList_t *lookup_res;
   int64_t hashval;
   int d1, d2;
   int found;
   hashval = link_hash_val(cur_link->end1_id, cur_link->end2_id, cur_link->compressedTails, hashtable->size);
   lookup_res = lookup_link(hashtable, cur_link->end1_id, cur_link->end2_id, cur_link->compressedTails ,cur_link->link_type, hashval);
   int cur_lib;
   
#ifdef DEBUG
   char linkAux1[LINK_CHUNK_SIZE];
   char linkAux2[LINK_CHUNK_SIZE];
   char tail2, tail1;

#endif
   
   if (lookup_res == NULL) {
      /* Should reserve a bucket from the appropriate heap and insert that to the table */
      if (myLog != 0) fprintf(myLog, "Thread %d: Storing %ld at hash[%ld]\n", MYTHREAD, *link_heap_pos, hashval);
      new_bucket = (linkList_t*) (linkList_heap + (*link_heap_pos));
      (*link_heap_pos)++;
      new_bucket->end1_id = cur_link->end1_id;
      new_bucket->end2_id = cur_link->end2_id;
      new_bucket->compressedTails = cur_link->compressedTails;
      new_bucket->link_type = cur_link->link_type;
      new_bucket->o1_length = cur_link->o1_length;
      new_bucket->o2_length = cur_link->o2_length;
      new_bucket->data = NULL;
      new_bucket->next = hashtable->table[hashval];
      hashtable->table[hashval] = new_bucket;

      lookup_res = new_bucket;
#ifdef DEBUG
      if (new_bucket->compressedTails == 0) {
         tail1 =3;
         tail2 =3;
      }
      if (new_bucket->compressedTails == 1) {
         tail1 =3;
         tail2 =5;
      }
      if (new_bucket->compressedTails == 2) {
         tail1 =5;
         tail2 =3;
      }
      if (new_bucket->compressedTails == 3) {
         tail1 =5;
         tail2 =5;
      }
      
      sprintf_chk(linkAux1, LINK_CHUNK_SIZE, "Contig%d.%d", new_bucket->end1_id, tail1);
      sprintf_chk(linkAux2, LINK_CHUNK_SIZE, "Contig%d.%d", new_bucket->end2_id, tail2);
      if (myLog != NULL) fprintf(myLog, "%s<=>%s\n", linkAux1, linkAux2);
#endif
      (*nLinks)++;
   } else {
      /* For SPANs check if alreads data with d1.d2 exist -- Equivalent to redundancy check */
      if ( cur_link->link_type == SPAN ) {
         span_data = lookup_res->data;
         d1 = cur_link->d1;
         d2 = cur_link->d2;
         while (span_data != NULL) {
            /* Do not add the same datum, i.e. same d1 and d2 since it is redundant */
            if ( (span_data->d1 == d1) && (span_data->d2 == d2) ) {
               cur_lib = cur_link->lib_id;
               bmas[cur_lib].libPairSummary[REDUN] += 1;
               (span_data->count)++;
               return 0;
            }
            span_data = span_data->next;
         }
      }
   }
   
   /* Add new datum to bucket's chain */
   if (myLog != NULL) fprintf(myLog, "Thread %d: adding new datum at %ld\n", MYTHREAD, *datum_heap_pos);
   new_datum = (datum_t*) (datum_heap + (*datum_heap_pos));
   (*datum_heap_pos)++;
   new_datum->next = lookup_res->data;
   lookup_res->data = new_datum;
   new_datum->endSeparation = cur_link->endSeparation;
   new_datum->lib_id = cur_link->lib_id;
   if (cur_link->link_type == SPAN) {
      if (myLog != NULL) fprintf(myLog, "Thread %d: adding SPAN link for lib %d\n", MYTHREAD, cur_link->lib_id);
      new_datum->d1 = cur_link->d1;
      new_datum->d2 = cur_link->d2;
      new_datum->count = 1;
      cur_lib = cur_link->lib_id;
      bmas[cur_lib].libPairSummary[ACCPT] += 1;
   }
   if (myLog != NULL) fprintf(myLog, "Thread %d: added_link_data\n", MYTHREAD);

   return 1;
}

/* Free routines */
int free_local_heaps(datum_t *datum_heap, linkList_t *linkList_heap)
{
   free(datum_heap);
   free(linkList_heap);
   
   return 0;
}

int free_hash_table(link_hash_table_t *hashtable)
{
   free(hashtable->table);
   
   return 0;
}
