#ifndef BUBBLE_FINDER_H
#define BUBBLE_FINDER_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <assert.h>
#include <upc.h>

#include "../../common/upc_compatibility.h"

#define BS 1

#ifndef KMER_LENGTH
#define KMER_LENGTH 51
#endif

#define PLUS '+'
#define MINUS '-'

#define USED_EXT 42
#define UNUSED_EXT 0


#define B_PLUS '0'
#define B_MINUS '1'
#define C_PLUS '2'
#define C_MINUS '3'

#define MAX_LINE_SIZE 1000
#define MAX_CONTIG_SIZE 900000

#define UNVISITED 0
#define VISITED 1

#define MAX_DEPTH 1

typedef struct stream_t stream_t;
struct stream_t{
   int objectId;
   char type;
   char strand;
};

/* Data_t data structure for tipInfo */
typedef struct data_t data_t;
struct data_t{
   int contig;
   char rcFlag;
   shared[] data_t *next;
};

/* Enhanced data_t data structure for tipInfo */
typedef struct enhanced_data_t enhanced_data_t;
struct enhanced_data_t{
   int contig;
   int length;
   char rcFlag;
   double depth;
   char *seq;
};

/* pdata_t data structure for pointsTo */
typedef struct pdata_t pdata_t;
struct pdata_t{
   int object;
   char type;
   shared[] pdata_t *next;
};

/* Hash-table entries data structure */
typedef struct list_t list_t;
struct list_t{
   char tipKey[2*KMER_LENGTH];
   int nContigs;
   shared[] list_t *next;
   shared[] data_t *data;
};

/* Hash-table entries data structure for points to */
typedef struct plist_t plist_t;
struct plist_t{
   char tipKey[KMER_LENGTH];
   int nContigs;
   shared[] plist_t *next;
   shared[] pdata_t *data;
};

/* Bucket data structure */
typedef struct bucket_t bucket_t;
struct bucket_t{
   upc_lock_t *bucket_lock;						// Lock to avoid race conditions during updates
   shared[] list_t *head;							// Pointer to the first entry of the hashtable
};

/* Bucket data structure for pointsTo */
typedef struct pbucket_t pbucket_t;
struct pbucket_t{
   upc_lock_t *bucket_lock;						// Lock to avoid race conditions during updates
   shared[] plist_t *head;							// Pointer to the first entry of the hashtable
};

/* Hash-table data structure */
typedef struct hash_table_t hash_table_t;
struct hash_table_t {
   int64_t size;										// Size of the hashtable
   shared[BS] bucket_t *table;					// Entries of the hashtable are pointers to buckets
};

/* Hash-table data structure for pointsTo */
typedef struct phash_table_t phash_table_t;
struct phash_table_t {
   int64_t size;										// Size of the hashtable
   shared[BS] pbucket_t *table;					// Entries of the hashtable are pointers to buckets
};

/* Memory heap data structure */
typedef struct memory_heap_t memory_heap_t;
struct memory_heap_t {
   shared[] list_t *my_heap_ptr;					// Pointer to my heap
   shared[] data_t *my_data_ptr;					// Pointer to my data
   int64_t my_heap_pos;
   int64_t my_data_pos;
};

/* Memory heap data structure for pointsTo */
typedef struct pmemory_heap_t pmemory_heap_t;
struct pmemory_heap_t {
   shared[] plist_t *my_heap_ptr;					// Pointer to my heap
   shared[] pdata_t *my_data_ptr;					// Pointer to my data
   int64_t my_heap_pos;
   int64_t my_data_pos;
};

/* cea data structure */
typedef struct cea_t cea_t;
struct cea_t {
   int contigID;
   int contigLen;
   char prevCodeL;
   char prevCodeR;
   char prevBase;
   char firstMer[KMER_LENGTH];
   char lastMer[KMER_LENGTH];
   char nextBase;
   char nextCodeL;
   char nextCodeR;
   char statusField;
   int nMers;
   double meanDepth;
};

/* extensions data structure */
typedef struct extensions_t extensions_t;
struct extensions_t {
   char data[2*KMER_LENGTH];
   char used_flag;
};

/* bubbletig data structure */
typedef struct bubbletig_t bubbletig_t;
struct bubbletig_t {
   char tipKey[2*KMER_LENGTH];
};

typedef struct sequence_t sequence_t;
struct sequence_t {
   shared[] char *data;
   int seqLength;
};

typedef struct local_sequence_t local_sequence_t;
struct local_sequence_t {
   char *data;
   int seqLength;
};

/* Split a merDepth line */
int splitDEPTH(char *inputLine, int *contigID, int *nMers, double *depth){
   char *token;
   char *aux;
   
   // merDepth lines are printed with the following command:
   // fprintf(myOutputFile, "Contig%d\t%d\t%f\n", contigID, nMers, mean);

   token = strtok_r(inputLine, "\t", &aux);
   (*contigID) = atoi(token+6);
   token = strtok_r(NULL, "\t", &aux);
   (*nMers) = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   (*depth) = atof(token);
   
   return 1;
}

/* Split a CEA line and store the contigs in a cea data structure */
int splitCEA(char *inputLine, cea_t *ceaEntry){
   char *token;
   char *aux;
   
   // CEA lines are printed with the following command:
   // fprintf(myOutputFile,"Contig%d\t[%c%c]\t(%c)\t%s\t%d\t%s\t(%c)\t[%c%c]\n", contigID, prevCodeL, prevCodeR, prevBase, firstKmer, contigLen, lastKmer, nextBase, nextCodeL, nextCodeR );
   
   token = strtok_r(inputLine, "\t", &aux);
   ceaEntry->contigID = atoi(token+6);
   token = strtok_r(NULL, "\t", &aux);
   ceaEntry->prevCodeL = *(token+1);
   ceaEntry->prevCodeR = *(token+2);
   token = strtok_r(NULL, "\t", &aux);
   ceaEntry->prevBase = *(token+1);
   token = strtok_r(NULL, "\t", &aux);
   memcpy(ceaEntry->firstMer, token, KMER_LENGTH*sizeof(char));
   token = strtok_r(NULL, "\t", &aux);
   ceaEntry->contigLen = atoi(token);
   token = strtok_r(NULL, "\t", &aux);
   memcpy(ceaEntry->lastMer, token, KMER_LENGTH*sizeof(char));
   token = strtok_r(NULL, "\t", &aux);
   ceaEntry->nextBase = *(token+1);
   token = strtok_r(NULL, "\t", &aux);
   ceaEntry->nextCodeL = *(token+1);
   ceaEntry->nextCodeR = *(token+2);
   
   return 1;
}

int isACGT(const char c) {
   if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
      return 1;
   else
      return 0;
}

char reverseComplementBaseStrict(const char base) {
   char rc;
   switch(base) {
      case 'A':
         rc = 'T'; break;
      case 'C':
         rc = 'G'; break;
      case 'G':
         rc = 'C'; break;
      case 'T':
         rc = 'A'; break;
      case '0':
         rc = '0'; break;
      default:
         fprintf(stderr, "unexpected base in revereseComplementBaseStrict: %c %d\n", base, (int) base);
         assert(0);
   }
   return rc;
}

void reverseComplementSeq(const char *seq, char *rc_seq, size_t seqLen)
{
   
   int end = seqLen-1;
   int start = 0;
   char temp;
   
   if(seqLen % 2 == 1) {
      /* has odd length, so handle middle */
      rc_seq[(start+end)/2] = reverseComplementBaseStrict( seq[(start+end)/2] );
   }
   while( start < end ) {
      
      temp = seq[end]; // in case this is inplace! (seq == rc_seq)
      rc_seq[end] = reverseComplementBaseStrict( seq[start] );
      rc_seq[start] = reverseComplementBaseStrict( temp );
      
      ++start;
      --end;
   }
}

void reverseComplementINPLACE(char *subcontig, int64_t size)
{
   reverseComplementSeq(subcontig, subcontig, size);
}

/* qsort struct comparison function (based on sequence length) */
int struct_cmp_by_length(const void *a, const void *b)
{
   enhanced_data_t *ia = (enhanced_data_t *)a;
   enhanced_data_t *ib = (enhanced_data_t *)b;
   return (int)(ia->length - ib->length);
}

int int64_comp(const void *a, const void *b)
{
   int64_t *ia = (int64_t *)a;
   int64_t *ib = (int64_t *)b;
   return (int)((*ia) - (*ib));
}

#endif //BUBBLE_FINDER_H

