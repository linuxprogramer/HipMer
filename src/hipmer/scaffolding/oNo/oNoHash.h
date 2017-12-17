#define SPLINT 0
#define SPAN 1
#define LINK_CHUNK_SIZE 100
#define BS 1
#define PLUS 0
#define MINUS 1
#include <upc.h>
#include <upc_nb.h>

char *myitoa2(int value, char* result, int base) {
   // check that the base if valid
   if (base < 2 || base > 36) { *result = '\0'; return result; }
   
   char* ptr = result, *ptr1 = result, tmp_char;
   int tmp_value;
   
   do {
      tmp_value = value;
      value /= base;
      *ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
   } while ( value );
   
   // Apply negative sign
   if (tmp_value < 0) *ptr++ = '-';
   *ptr-- = '\0';
   while(ptr1 < ptr) {
      tmp_char = *ptr;
      *ptr--= *ptr1;
      *ptr1++ = tmp_char;
   }
   return result;
}

typedef struct oNolink_t oNolink_t;
struct oNolink_t{
   int end1_id;                  // ID of link's endpoint 1
   int end2_id;                  // ID of link's endpoint 2
   char end1_s;
   char end2_s;
   char link_type;      // Type of link: SPLINT or SPAN
   char lib_id;                  // Library id
   int f1;
   int f2;
   int f3;
   int f4;
   int f5;
};

typedef shared[] oNolink_t* oNoShared_link_ptr;

/* Link heap data structure */
typedef struct oNoLink_heap_t oNoLink_heap_t;
struct oNoLink_heap_t {
	shared[BS] oNoShared_link_ptr *link_ptr;				// Pointers to shared memory heaps
	shared[BS] int64_t *heap_indices;						// Indices of remote heaps
};

/* Creates shared heaps that required for links hash table - IMPORTANT: This is a collective function */
int create_oNolink_heaps(int64_t size, oNoLink_heap_t *link_heap)
{
	link_heap->link_ptr = (shared[BS] oNoShared_link_ptr*) upc_all_alloc(THREADS, sizeof(oNoShared_link_ptr));
   link_heap->heap_indices = (shared[BS] int64_t*) upc_all_alloc(THREADS, sizeof(int64_t));
	link_heap->link_ptr[MYTHREAD] = (shared[] oNolink_t*) upc_alloc( (size+THREADS)/THREADS * sizeof(oNolink_t));
	link_heap->heap_indices[MYTHREAD] = 0;
	upc_barrier;
	
	return 0;
}

/* Allocate local arrays used for book-keeping when using aggregated upc_memputs */
int allocate_oNolocal_buffs(oNolink_t **local_buffs, int64_t **local_index)
{
	(*local_buffs) = (oNolink_t*) malloc_chk(THREADS * LINK_CHUNK_SIZE * sizeof(oNolink_t));
	(*local_index) = (int64_t*) malloc_chk(THREADS * sizeof(int64_t));
	memset((*local_index), 0, THREADS * sizeof(int64_t));

   return 0;
}

/* Adds a link to the shared link heap */
int add_oNolink_to_shared_heaps(oNolink_t *new_entry, int64_t hashval, int64_t *local_index, oNolink_t *local_buffs, oNoLink_heap_t link_heap)
{
   int64_t store_pos;
   int remote_thread = hashval % THREADS;
   
   /* Store link first to local buffer designated for remote thread */
   if (local_index[remote_thread] <= LINK_CHUNK_SIZE - 1) {
      memcpy( &(local_buffs[local_index[remote_thread] + remote_thread * LINK_CHUNK_SIZE]), new_entry, sizeof(oNolink_t));
      local_index[remote_thread]++;
   }
   
   /* If buffer for that thread is full, do a remote upc_memput() */
   if (local_index[remote_thread] == LINK_CHUNK_SIZE) {
      store_pos = bupc_atomicI64_fetchadd_strict(&(link_heap.heap_indices[remote_thread]), LINK_CHUNK_SIZE);
      upc_memput( (shared[] oNolink_t*)  ((link_heap.link_ptr[remote_thread]) + store_pos)  , &(local_buffs[remote_thread * LINK_CHUNK_SIZE]), (LINK_CHUNK_SIZE) * sizeof(oNolink_t));
      local_index[remote_thread] = 0;
   }

   return 0;
}

/* Adds remaining links to the shared heaps. Should be called when all calls "add_oNolink_to_shared_heaps()" have been done */
int add_rest_oNolinks_to_shared_heaps(int64_t *local_index, oNolink_t *local_buffs, oNoLink_heap_t link_heap)
{
   int i;
   int64_t store_pos;
   
   for (i = 0; i < THREADS; i++) {
		if (local_index[i] != 0) {
			store_pos = bupc_atomicI64_fetchadd_strict(&(link_heap.heap_indices[i]), local_index[i]);
			upc_memput( (shared[] oNolink_t*)  ((link_heap.link_ptr[i]) + store_pos)  , &(local_buffs[i * LINK_CHUNK_SIZE]), (local_index[i]) * sizeof(oNolink_t));
		}
	}
   
   return 0;
}

/*
int64_t oNolinkhash(int64_t linkhash_size, char *linkname)
{
   unsigned long hashval;
   hashval = 5381;
   for(; *linkname != '\0'; linkname++) hashval = (*linkname) +  (hashval << 5) + hashval;
   return (hashval % linkhash_size);
}
*/

/* Data structure to store data regarding a link */
typedef struct oNoDatum_t oNoDatum_t;
struct oNoDatum_t{
   char lib_id;
   int f1;
   int f2;
   int f3;
   int f4;
   int f5;
   char link_type;      // Type of link: SPLINT or SPAN
	oNoDatum_t *next;                // Pointer to next entry in the same datum list
};

/* Essentially the key now is the combination: end1_id | end2_id | linktype */
typedef struct oNolinkList_t oNolinkList_t;
struct oNolinkList_t{
   int end1_id;                  // ID of link's endpoint 1
   int end2_id;                  // ID of link's endpoint 2
   char end1_s;
   char end2_s;
   oNoDatum_t *data;                // Data list
	oNolinkList_t *next;					// Pointer to next entry in the same bucket
};

typedef struct oNolink_hash_table_t_ {
   int64_t size;                     /* the size of the table */
   oNolinkList_t **table;           /* the table elements */
} oNolink_hash_table_t;

/* Create a link_hash_table */
oNolink_hash_table_t *create_oNolink_hash_table(int64_t size)
{
   oNolink_hash_table_t *new_table;
   int64_t i;
   
   if (size<1) return NULL;
   
   /* Attempt to allocate memory for the table structure */
   if ((new_table = (oNolink_hash_table_t*) malloc_chk(sizeof(oNolink_hash_table_t))) == NULL) {
      return NULL;
   }
   
   /* Attempt to allocate memory for the table itself */
   if ((new_table->table = (oNolinkList_t **) malloc_chk(sizeof(oNolinkList_t*) * size)) == NULL) {
      return NULL;
   }
   
   /* Initialize the elements of the table */
   for(i=0; i<size; i++) new_table->table[i] = NULL;
   
   /* Set the table's size */
   new_table->size = size;
   
   return new_table;
}

int64_t oNolink_hash_val(int e1, int e2, char s1, char s2 , int64_t linkhash_size)
{
   char linkname[100];
   char s1c;
   char s2c;
   
   //s1c = (s1 == 3) ? '3' : '5';
   //s2c = (s2 == 3) ? '3' : '5';
   
   if (s1 == 3) {
      s1c = '3';
   } else {
      s1c = '5';
   }
   
   if (s2 == 3) {
      s2c = '3';
   } else {
      s2c = '5';
   }
   
   sprintf(linkname, "%d.%c<=>%d.%c" ,e1,s1c,e2,s2c);
   
   return hashstr(linkhash_size, linkname);
}

/* Lookup function for link hash table */
oNolinkList_t *lookup_oNolink(oNolink_hash_table_t *hashtable, int e1, int e2, char s1 ,char s2, int64_t hashval)
{
   oNolinkList_t *list;
   
   for(list = hashtable->table[hashval]; list != NULL; list = list->next) {
      if ((e1 == list->end1_id) && (e2 == list->end2_id) && (s1 == list->end1_s) && (s2 == list->end2_s)) return list;
   }

   return NULL;
}

/* Implement efficient memory allocator to avoid multiple malloc() calls */
int create_oNolocal_heaps(int64_t max_size, oNoDatum_t **datum_heap, oNolinkList_t **linkList_heap)
{
   (*datum_heap) = (oNoDatum_t*) malloc_chk(max_size * sizeof(oNoDatum_t));
   (*linkList_heap) = (oNolinkList_t*) malloc_chk(max_size * sizeof(oNolinkList_t));
   
   return 0;
}

/* Add link data to the hashtable */
int add_oNolink_data(oNolink_hash_table_t *hashtable, oNolink_t *cur_link, int64_t *datum_heap_pos, oNoDatum_t *datum_heap, int64_t *link_heap_pos, oNolinkList_t *linkList_heap)
{
   oNolinkList_t *new_bucket;
   oNoDatum_t *new_datum, *span_data;
   oNolinkList_t *lookup_res;
   int64_t hashval;
   int d1, d2;
   hashval = oNolink_hash_val(cur_link->end1_id, cur_link->end2_id, cur_link->end1_s, cur_link->end2_s, hashtable->size);
   lookup_res = lookup_oNolink(hashtable, cur_link->end1_id, cur_link->end2_id, cur_link->end1_s, cur_link->end2_s, hashval);
   
   if (lookup_res == NULL) {
      /* Should reserve a bucket from the appropriate heap and insert that to the table */
      new_bucket = (oNolinkList_t*) (&linkList_heap[(*link_heap_pos)]);
      (*link_heap_pos)++;
      new_bucket->end1_id = cur_link->end1_id;
      new_bucket->end2_id = cur_link->end2_id;
      new_bucket->end1_s = cur_link->end1_s;
      new_bucket->end2_s = cur_link->end2_s;
      //new_bucket->link_type = cur_link->link_type;
      new_bucket->data = NULL;
      new_bucket->next = hashtable->table[hashval];
      hashtable->table[hashval] = new_bucket;
      lookup_res = new_bucket;
   }

   
   /* Add new datum to bucket's chain */
   new_datum = (oNoDatum_t*) (datum_heap + (*datum_heap_pos));
   (*datum_heap_pos)++;
   new_datum->next = lookup_res->data;
   lookup_res->data = new_datum;
   new_datum->lib_id = cur_link->lib_id;
   new_datum->f1 = cur_link->f1;
   new_datum->f2 = cur_link->f2;
   new_datum->f3 = cur_link->f3;
   new_datum->link_type = cur_link->link_type;
   if (cur_link->link_type == SPAN) {
      new_datum->f5 = cur_link->f5;
   } else if (cur_link->link_type == SPLINT) {
      new_datum->f4 = cur_link->f4;
   }
   
   return 0;
}

/* Free routines */
int free_oNolocal_heaps(oNoDatum_t *datum_heap, oNolinkList_t *linkList_heap)
{
   free(datum_heap);
   free(linkList_heap);
   
   return 0;
}

int free_oNohash_table(oNolink_hash_table_t *hashtable)
{
   free(hashtable->table);
   
   return 0;
}
