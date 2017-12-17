#include <upc.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>

#define TRUE 1
#define FALSE 0
#define UU_CONTIG 0
#define NON_UU_CONTIG 1
#define BS 1

/* Conitg data structure */
typedef struct contig_t contig_t;
struct contig_t{
	shared[] char *contig;						// Actual content of the contig
	int64_t length;
   int64_t contig_id;
   int64_t uuType;
   int64_t parentContigID;
   int64_t offsetInParent;
   int64_t parentLength;
};

typedef shared[] contig_t* contig_t_ptr;

/* Data structure to hold a contig pointer and posInContig */
typedef struct contigInfo contigInfo;
struct contigInfo{
   shared[] contig_t *my_contig;
   int32_t posInContig;
};

/* Hash-table entries, i.e. seeds */
typedef shared[] contigInfo * contigInfoArrayPtr_t;
typedef struct list_t list_t;
struct list_t{
   contigInfoArrayPtr_t contigInfoArrayPtr;        // Pointer to the list of contigInfo entries
   shared[] list_t *next;                          // Pointer to next entry in the same bucket
#ifdef STORE_OPT
   int16_t nExtraContigs;
   int16_t posInSubArray;
#else
   int64_t nExtraContigs;                              // Number of extra contigs this seed can be found
   int64_t posInSubArray;                          // last written record in contigInfoArrayPtr
#endif
   contigInfo firstContig;                         // Information for the first contig
   unsigned char packed_key[KMER_PACKED_LENGTH];	// The packed key of the seed
};

typedef shared[] list_t* shared_heap_ptr;

/* Bucket data structure */
typedef struct bucket_t bucket_t;
struct bucket_t{
#ifndef STORE_OPT
        upc_lock_t *bucket_lock;
#endif
	shared[] list_t *head;							// Pointer to the first entry of the hashtable
};

/* Hash-table data structure */
typedef struct hash_table_t hash_table_t;
struct hash_table_t {
	int64_t size;										// Size of the hashtable
	shared[BS] bucket_t *table;					// Entries of the hashtable are pointers to buckets
};

/* Memory heap data structure */
typedef struct memory_heap_t memory_heap_t;
struct memory_heap_t {
	shared[BS] shared_heap_ptr *heap_ptr;		       // Pointers to shared memory heaps
	shared[BS] int64_t *heap_indices;		       // Indices of remote heaps
        shared[BS] int64_t *heap_ptr_sizes;                    // maximum allocated length of heap_ptr(s)
        shared[BS] contigInfoArrayPtr_t *contigInfoArrays;     // allocation for extraContigs
};

