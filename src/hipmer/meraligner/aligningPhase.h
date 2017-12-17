#ifndef __ALIGNING_PHASE_H
#define __ALIGNING_PHASE_H

#define PLUS 0
#define MINUS 1
#define ON 1
#define OFF 0
#define SUCCESS 1
#define FAIL 0
#define MAX_READSIZE 102
#define ALIGNS_SLACK 10

#define MAX_ALIGNMENT_LINE_SIZE 200

#include <assert.h>
#include <upc_collective.h>
#include "../common/upc_compatibility.h"
#include "ssw.h"
#include "cache_infrastructure.h"
#include "../fqreader/fq_reader.h"
#include "../common/common.h"
#include "../common/Buffer.h"

#ifdef PROFILE
extern double readmap_time;
extern double fetch_contigs_time;
extern double hash_table_lookups;
extern double input_IO;
extern double output_IO;
#endif

int read_map(FILE *resultFd, hash_table_t *hashtable, Buffer cur_read_buf, Buffer read_info_buf, Buffer quals, shared int64_t *cachePtr, shared[1] contigDataPtr *cacheTableContig, int min_contig_length, int8_t *nt_table, int8_t *mat, FILE *logFD, int chunk_size, int64_t cacheCapacity, shared[1] list_t *cacheTable, shared[1] int64_t *cacheFlags);


static shared int _all_max_read_len[THREADS];

int64_t parallelAligner(hash_table_t *hashtable, char *filename, FILE *res_fd1, FILE *res_fd2, int filesPerPair, int readNumber, double estimatedInsertSize, double estimatedSigma, int totalContigs, int64_t cacheSizeInMB, int min_contig_length, int chunk_size, int64_t cacheCapacity, shared[1] list_t *cacheTable, shared[1] int64_t *cacheFlags, shared int64_t *cachePtr, shared[1] contigDataPtr *cacheTableContig, int64_t *readsMapped, int64_t *readsProcessed, const char *base_dir, const char *libName)
{
   long total_read;
   int idlength, read_length;
   int j, k, found_flag, insert_size, err_code, pos1, pos2, orientation, test1, test2, test3, strand1, strand2, sStart1, sStart2, sStop1, sStop2, sLength1, min, max, extent, safeDistance, leftBoundary, rightBoundary;
   int64_t i, partial_res, readsRead, l, m;
   int locs[4];
   int64_t alignmentsFound = 0, curReadAlignments = 0;
   UPC_TICK_T start_map, end_map;
   double map_time = 0.0;
   double start_readmap_timer, end_readmap_timer;
   double start_timer, end_timer;
   int finished;
   int match = 1, mismatch = 3, gap_open = 5, gap_extension = 2, ambiguity_score = 2;
   int whichOutput = 0;
   FILE *logFD = NULL;
#ifdef DEBUG
   char logfile_name[MAX_FILE_PATH];
   sprintf(logfile_name,"parallelAligner_logfile_%d", MYTHREAD);
   get_rank_path(logfile_name, MYTHREAD);
   logFD = fopen_chk(logfile_name,"a");
   fprintf(stderr, "Thread %d: Opening %s for more detailed logs\n", MYTHREAD, logfile_name);
   fprintf(logFD, "Starting parallelAligner\n");
#endif
   
   /* This table is used to transform nucleotide letters into numbers. */
   int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};
   
   // initialize scoring matrix for genome sequences
	//  A  C  G  T	N (or other ambiguous code)
	//  2 -2 -2 -2 	0	A
	// -2  2 -2 -2 	0	C
	// -2 -2  2 -2 	0	G
	// -2 -2 -2  2 	0	T
	//	0  0  0  0  0	N (or other ambiguous code)
	int8_t* mat = (int8_t*) calloc_chk(25, sizeof(int8_t));
	for (l = k = 0; l < 4; ++l) {
		for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = -ambiguity_score; // ambiguous base: no penalty
	}
	for (m = 0; m < 5; ++m) mat[k++] = -ambiguity_score;
   
   upc_barrier;
   if (filesPerPair == 2) {
     whichOutput = readNumber;
   }
   
   total_read = 0;

   int from_shm = 0;
   if (strstr(base_dir, "/dev/shm"))
       from_shm = 1;

   fq_reader_t fqr = create_fq_reader();
   open_fq(fqr, filename, from_shm);

   if (MYTHREAD == 0) printf("Start with file %s\n", filename);

   Buffer id = initBuffer(MAX_READ_NAME_LEN);
   Buffer seq = initBuffer(MAX_READ_LEN), quals = initBuffer(MAX_READ_LEN);

#ifdef DEBUG
   fprintf(logFD, "Thread %d: Starting to read %s\n", MYTHREAD, filename);
#endif

   while (1) {

#ifdef PROFILE
       start_timer = UPC_TICKS_NOW();
#endif
       int found_next = get_next_fq_record(fqr, id, seq, quals);
#ifdef PROFILE
       end_timer = UPC_TICKS_NOW();
       input_IO += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
       if (!found_next)
           break;
       total_read++;

#ifdef PROFILE
       start_readmap_timer = UPC_TICKS_NOW();
#endif
      
#ifdef DEBUG
       fprintf(logFD, "Thread %d: processing %s\n", MYTHREAD, getStartBuffer(id));
#endif
       
       curReadAlignments = read_map(whichOutput == 0 ? res_fd1 : res_fd2, hashtable, seq, id, quals, cachePtr, cacheTableContig, min_contig_length, nt_table, mat, logFD, chunk_size, cacheCapacity, cacheTable, cacheFlags);
       if (filesPerPair == 1) {
           whichOutput = (whichOutput+1) % 2;
       }
       alignmentsFound += curReadAlignments;
       if (curReadAlignments > 0) {
           (*readsMapped)++;
       }
         
#ifdef PROFILE
       end_readmap_timer = UPC_TICKS_NOW();
       readmap_time += (UPC_TICKS_TO_SECS(end_readmap_timer-start_readmap_timer));
#endif
   }

   close_fq(fqr);


   // reduce to max read length over all threads and write to file
   shared int *per_thread_max_read_len = upc_all_alloc(THREADS, sizeof(int));
   per_thread_max_read_len[MYTHREAD] = fqr->max_read_len;
   bupc_all_reduce_allI(_all_max_read_len, per_thread_max_read_len, UPC_MAX, THREADS, 1, NULL, 
                        UPC_IN_MYSYNC|UPC_OUT_MYSYNC);
   upc_all_free(per_thread_max_read_len);

   if (!MYTHREAD) {
      printf("Max read length for lib %s: %d\n", libName, _all_max_read_len[0]);
      char buf[MAX_FILE_PATH];
      sprintf(buf, "%s-readlen.txt", libName);
      printf("about to write to %s\n", buf);
      FILE *f = fopen_rank_path(buf, "a", -1);
      fprintf(f, "%d\n", _all_max_read_len[0]);
      fclose(f);
   }

   
   (*readsProcessed) += total_read;
   freeBuffer(id); id = NULL;
   freeBuffer(seq); seq = NULL;
   freeBuffer(quals); quals = NULL;
   

#ifdef DEBUG
   fprintf(logFD, "Thread %d: Done processing: %ld\n", MYTHREAD, total_read);
   fclose(logFD);
#endif
   
   upc_barrier;
   if (MYTHREAD == 0) {
      printf("Done with file %s, total reads read are %ld\n", filename, total_read);
   }
   return alignmentsFound;
}

int read_map(FILE *resultFd, hash_table_t *hashtable, Buffer cur_read_buf, Buffer read_info_buf, Buffer quals, shared int64_t *cachePtr, shared[1] contigDataPtr *cacheTableContig, int min_contig_length, int8_t *nt_table, int8_t *mat, FILE *logFD, int chunk_size, int64_t cacheCapacity, shared[1] list_t *cacheTable, shared[1] int64_t *cacheFlags)
{
   int front_ptr, prev_front_contig_id, prev_front_strand, front_contig_id, back_contig_id, new_front_contig_id, back_ptr, matched_pos_cont, matched_pos_read, commun_opt_flag, m, k;
   int lex_ind, contig_length, found_flag, bypass_opt, strand, contig_fw, has_been_flipped = OFF;
   int forward_search, interest_zone;
   int posInContig, posInContigPrimer, effective_length;
   int offsetssw;
   shared[] list_t *lookup_res;
   list_t seedEntry;
   int foundSeed;
   shared[] contig_t *cur_contig;
   contig_t local_contig;
   char *kmer_to_search;
   const char *plus_strand="Plus", *minus_strand="Minus";
   const char *print_strand;
   char cur_kmer[KMER_LENGTH+1];
   char cur_kmer2[KMER_LENGTH+1];
   char rc_cur_kmer[KMER_LENGTH+1];
   rc_cur_kmer[KMER_LENGTH] = '\0';
   cur_kmer[KMER_LENGTH] = '\0';
   cur_kmer2[KMER_LENGTH] = '\0';
   char *search_ptr;
   int contigs_matched = 0, min_len;

   // quals is not used, so repurpose for contig_copy
   Buffer contig_copy_buf = quals;
   resetBuffer(contig_copy_buf);
   char *contig_copy = NULL;

   double start_fetching_timer, end_fetching_timer;
   double start_timer, end_timer;
   int8_t *num = NULL, *ref_num = NULL;
   const char *read_seq, *ref_seq;
   s_align *r;
   s_profile *profile = NULL;
   int match = 1, mismatch = 3, gap_open = 5, gap_extension = 2;
   int *read_alignments = NULL, *read_alignments_scores = NULL;
   char *alignment_buffer = NULL;
   int ind;
   int bound, found;
   int contigType;

   contigDataPtr remote_contig_ptr;
   
   /* Copy input read and input read identifier in valid string buffers */
   const char *cur_read = getStartBuffer(cur_read_buf);
   int read_length = getLengthBuffer(cur_read_buf);
   assert(read_length > 0);
   
   const char *read_info = getStartBuffer(read_info_buf);
   int read_info_length = getLengthBuffer(read_info_buf);
   assert(read_info_length > 0);
   
   /* Lookup for the first "valid" seed in the read */
   front_ptr = 0;
   lookup_res = NULL;
   foundSeed = 0;
   while ( (lookup_res == NULL) && (front_ptr+KMER_LENGTH <= read_length ) ) {
      memcpy(cur_kmer, cur_read+front_ptr, KMER_LENGTH * sizeof(char));
      if (strchr(cur_kmer,'N') != NULL) {
         front_ptr++;
      } else {
		 reverseComplementKmer(cur_kmer, rc_cur_kmer);
         lex_ind = strcmp(cur_kmer, rc_cur_kmer);
         kmer_to_search = ( lex_ind < 0) ? cur_kmer : rc_cur_kmer;
#ifdef PROFILE
         start_timer = UPC_TICKS_NOW();
#endif
         
#ifdef USE_SWCACHE
         lookup_res = lookupKmerInCache(hashtable, (unsigned char*) kmer_to_search, cacheCapacity, cacheTable, cacheFlags);
#else
         lookup_res = lookup_kmer(hashtable, (unsigned char*) kmer_to_search);
#endif
         
#ifdef PROFILE
         end_timer = UPC_TICKS_NOW();
         hash_table_lookups += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
         if (lookup_res == NULL) {
            /* The k-mer is not a valid seed */
            front_ptr++;
         } else {
            /* The k-mer has been added to a contig */
#ifdef PROFILE
            start_timer = UPC_TICKS_NOW();
#endif
            seedEntry = (list_t) *lookup_res;
            foundSeed = 1;
#ifdef PROFILE
            end_timer = UPC_TICKS_NOW();
            hash_table_lookups += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
         }
      }
   }
   
   if (foundSeed == 0) {
      //printf("The read contains no valid kmer - forward search\n");
      /* FIXME: Remove this "guard alignment" when processing a single read file containing both lanes /1 and /2 */
      /* Print foo alignment in result so that we know that no valid alignment is found */
      fprintf(resultFd, "MERALIGNER-F\t%s\t\t0\t0\tContig0\t0\t0\t0\tPlus\t0\t0\t0\t0\n", read_info);
#ifdef DEBUG
      fprintf(logFD, "No seed found for %s\n", read_info);
#endif
      
      return 0;
   }
   
   /* Copy remote contig to manipulate locally */
   cur_contig = seedEntry.firstContig.my_contig;
   upc_memget(&local_contig, cur_contig, sizeof(contig_t));
   contig_length = local_contig.length;
   
   front_contig_id = local_contig.contig_id;

#ifdef PROFILE
   start_fetching_timer = UPC_TICKS_NOW();
#endif
   
#ifdef USE_SWCACHE
   remote_contig_ptr = findRemoteContigPtr(front_contig_id, contig_length, cur_contig, 0, cachePtr, cacheTableContig);
#else
   remote_contig_ptr = (contigDataPtr) (&(cur_contig->contig[0]));
#endif
#ifdef DEBUG2
   fprintf(logFD, "Thread %d: Copying %d seq for contig %d %ld\n", MYTHREAD, contig_length, front_contig_id, local_contig.parentContigID);
#endif
   contig_copy = resetRawBuffer(contig_copy_buf, contig_length);
   upc_memget(contig_copy, remote_contig_ptr, contig_length * sizeof(char));
   assert(contig_copy[contig_length] == '\0');
   
#ifdef PROFILE
   end_fetching_timer = UPC_TICKS_NOW();
   fetch_contigs_time += (UPC_TICKS_TO_SECS(end_fetching_timer-start_fetching_timer));
#endif
   posInContig = seedEntry.firstContig.posInContig;
   contigType = local_contig.uuType;
   assert(posInContig >= 0 && contig_length >= posInContig + KMER_LENGTH);
   
   /* Check if we should align with respect to the current contig or the reverse complement */
   strand = PLUS;
   assert( getLengthBuffer(contig_copy_buf) >= posInContig + KMER_LENGTH );
   if ( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) != 0 ) {
      strand =  MINUS;
#ifdef DEBUG
      //fprintf(stderr, "Thread %d: revcomp3 contig %d, %d (%d len): %.*s\n", MYTHREAD, local_contig.parentContigID, front_contig_id, contig_length, contig_length, contig_copy);
#endif
      reverseComplementINPLACE(contig_copy, contig_length);
      posInContig = contig_length - (KMER_LENGTH + posInContig);
      assert( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) == 0);
   }
   
   /* Check if sufficient bases in contig such that it covers the whole read and if we have UU-contig */
   if ( ((contig_length - posInContig) >= (read_length - front_ptr)) && (front_ptr == 0) && (contigType == UU_CONTIG)) {
      if ( memcmp(cur_read + front_ptr, contig_copy + posInContig, (read_length - front_ptr) * sizeof(char)) == 0 ) {
         
         /* Translate coordinates to parent contig */
         contig_length = local_contig.parentLength;
         if (strand == PLUS) {
            posInContig = posInContig + local_contig.offsetInParent;
         } else {
            posInContig = posInContig + contig_length - (local_contig.length + local_contig.offsetInParent);
         }
         front_contig_id = local_contig.parentContigID;

         
         /* Bypass optimization check has been successful */
         print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
#ifdef PROFILE
         start_timer = UPC_TICKS_NOW();
#endif
         if (strand == PLUS) {
            fprintf(resultFd, "MERALIGNER-0\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, 1, read_length, read_length, front_contig_id, posInContig+1, posInContig+read_length, contig_length, print_strand );
         } else {
            fprintf(resultFd, "MERALIGNER-0\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, 1, read_length, read_length, front_contig_id, contig_length-(posInContig+read_length-1), contig_length - posInContig, contig_length, print_strand );
         }
#ifdef PROFILE
         end_timer = UPC_TICKS_NOW();
         output_IO += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
#ifdef DEBUG
         fprintf(logFD, "Found UU %s for %s\n", strand ? "plus" : "minus", read_info);
#endif
         return 1;
      }
   }
   
   /* Allocate alignment history buffers */
   int maxAlignments = (read_length - KMER_LENGTH + 1) * ALIGNS_SLACK;
   read_alignments = (int*) malloc_chk(maxAlignments * sizeof(int));
   read_alignments_scores = (int*) malloc_chk(maxAlignments * sizeof(int));
   alignment_buffer = (char*) malloc_chk(maxAlignments * MAX_ALIGNMENT_LINE_SIZE * sizeof(char));

   for (ind = 0; ind < maxAlignments; ind++) {
      read_alignments[ind] = -1;
      read_alignments_scores[ind] = 0;
   }
   bound = 1;

   /* Use SSW to align as much as we can from this read (align from start of the read to allow misalignments at the beginning)*/
   num = (int8_t*) malloc_chk(read_length * sizeof(int8_t));
   ref_num = (int8_t*) malloc_chk(contig_length * sizeof(int8_t));
   read_seq = cur_read;
   for (m = 0; m < read_length; ++m) num[m] = (int8_t) nt_table[(int)read_seq[m]];
   assert(num);
   profile = ssw_init(num, read_length, mat, 5, 2);
   
   /* Slack here should be equal to front_ptr  */
   posInContigPrimer = ( (posInContig - front_ptr - read_length) < 0 ) ? 0 : posInContig - front_ptr - read_length;
   ref_seq = contig_copy + posInContigPrimer;
   effective_length = contig_length - posInContigPrimer;
   for (m = 0; m < effective_length; ++m) ref_num[m] = (int8_t) nt_table[(int)ref_seq[m]];
   min_len = (read_length/2 >= 15) ? read_length / 2 : 15;
   assert(profile);
   r = ssw_align (profile, ref_num, effective_length, gap_open, gap_extension, 8, 0, 0, min_len);
   
   assert(ref_num);
   free(ref_num); ref_num = NULL;
   
   print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
#ifdef PROFILE
   start_timer = UPC_TICKS_NOW();
#endif
   
   /* Translate coordinates to parent contig */
   contig_length = local_contig.parentLength;
   if (strand == PLUS) {
      posInContigPrimer = posInContigPrimer + local_contig.offsetInParent;
   } else {
      posInContigPrimer = posInContigPrimer + contig_length - (local_contig.length + local_contig.offsetInParent);
   }
   front_contig_id = local_contig.parentContigID;
   
   read_alignments[0] = front_contig_id;
   read_alignments_scores[0] = r->score1;
   
   
   if (strand == PLUS) {
      sprintf(alignment_buffer, "MERALIGNER-1\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, r->ref_begin1+1+posInContigPrimer,  r->ref_end1+1+posInContigPrimer, contig_length, print_strand );
   } else {
      sprintf(alignment_buffer, "MERALIGNER-1\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, contig_length-(r->ref_end1+posInContigPrimer), contig_length - (r->ref_begin1+posInContigPrimer), contig_length, print_strand );
   }
#ifdef PROFILE
   end_timer = UPC_TICKS_NOW();
   output_IO += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
   
   matched_pos_read = r->read_end1;
   align_destroy(r);
   
   contigInfo *localContigInfoArray = NULL;
   int z;
   
   /* If we have multiple contigs that this seed is extracted from, find alignments with those as well */
   if ( seedEntry.nExtraContigs > 0 ) {
      localContigInfoArray = (contigInfo*) malloc_chk((seedEntry.nExtraContigs) * sizeof(contigInfo));
      assert(seedEntry.contigInfoArrayPtr != NULL);
#ifdef DEBUG
      fprintf(logFD, "Thread %d: getting %lld contigInfoArrayPointers for %s from %d\n", MYTHREAD, (long long) seedEntry.nExtraContigs, cur_kmer, (int) upc_threadof(seedEntry.contigInfoArrayPtr));
#endif
      assert( IS_VALID_UPC_PTR(seedEntry.contigInfoArrayPtr ) );
      upc_memget(localContigInfoArray, seedEntry.contigInfoArrayPtr , (seedEntry.nExtraContigs) * sizeof(contigInfo));
      for (z = 0; z < seedEntry.nExtraContigs; z++) {
         /* Copy remote contig to manipulate locally */
         cur_contig = localContigInfoArray[z].my_contig;
         upc_memget(&local_contig, cur_contig, sizeof(contig_t));
         contig_length = local_contig.length;
         contig_copy = resetRawBuffer(contig_copy_buf, contig_length);
         assert(contig_copy[contig_length] == '\0');
         front_contig_id = local_contig.contig_id;
         
#ifdef PROFILE
         start_fetching_timer = UPC_TICKS_NOW();
#endif
         
#ifdef USE_SWCACHE
         remote_contig_ptr = findRemoteContigPtr(front_contig_id, contig_length, cur_contig, 0, cachePtr, cacheTableContig);
#else
         remote_contig_ptr = (contigDataPtr) (&(cur_contig->contig[0]));
#endif
#ifdef DEBUG2
         fprintf(logFD, "Thread %d: Copying %d seq for contig %d %ld\n", MYTHREAD, contig_length, front_contig_id, local_contig.parentContigID);
#endif
         upc_memget(contig_copy, remote_contig_ptr, contig_length * sizeof(char));
         assert(contig_copy[contig_length] == '\0');
         
#ifdef PROFILE
         end_fetching_timer = UPC_TICKS_NOW();
         fetch_contigs_time += (UPC_TICKS_TO_SECS(end_fetching_timer-start_fetching_timer));
#endif
         posInContig = localContigInfoArray[z].posInContig;
         assert(posInContig >= 0 && contig_length >= posInContig + KMER_LENGTH);
         
         /* Check if we should align with respect to the current contig or the reverse complement */
         strand = PLUS;
         assert( getLengthBuffer(contig_copy_buf) >= posInContig + KMER_LENGTH ); 
         if ( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) != 0 ) {
            strand =  MINUS;
#ifdef DEBUG
         fprintf(logFD, "Thread %d: revcomp2 contig %d (%d len): %.*s\n", MYTHREAD, front_contig_id, contig_length, contig_length, contig_copy);
#endif
            reverseComplementINPLACE(contig_copy, contig_length);
            posInContig = contig_length - (KMER_LENGTH + posInContig);
            assert( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) == 0);
         }
         
         /* Use SSW to align as much as we can from this read (align from start of the read to allow misalignments at the beginning)*/
         assert(!ref_num);
         ref_num = (int8_t*) malloc_chk(contig_length * sizeof(int8_t));
         read_seq = cur_read;
         posInContigPrimer = ( (posInContig - front_ptr - read_length) < 0 ) ? 0 : posInContig - front_ptr - read_length;
         assert( posInContigPrimer < getLengthBuffer(contig_copy_buf) );
         ref_seq = contig_copy + posInContigPrimer;
         effective_length = contig_length - posInContigPrimer;
         for (m = 0; m < effective_length; ++m) ref_num[m] = (int8_t) nt_table[(int)ref_seq[m]];
         min_len = (read_length/2 >= 15) ? read_length / 2 : 15;
         assert(profile);
         assert(ref_num);
         r = ssw_align (profile, ref_num, effective_length, gap_open, gap_extension, 8, 0, 0, min_len);
         
         free(ref_num); ref_num = NULL;
         found = 0;
         front_contig_id = local_contig.parentContigID;

         for (ind = 0; ind < bound; ind++) {
            assert(ind < maxAlignments);
            /* Search if alignment with the same contig has been already found and check scores */
            if (read_alignments[ind] == front_contig_id) {
               /* Already have found an alignment with the same contig, check if the score is higer now and if so update output result */
               if (r->score1 > read_alignments_scores[ind]) {
                  read_alignments_scores[ind] = r->score1;
                  print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
                  
                  /* Translate coordinates to parent contig */
                  contig_length = local_contig.parentLength;
                  if (strand == PLUS) {
                     posInContigPrimer = posInContigPrimer + local_contig.offsetInParent;
                  } else {
                     posInContigPrimer = posInContigPrimer + contig_length - (local_contig.length + local_contig.offsetInParent);
                  }
                  front_contig_id = local_contig.parentContigID;
                  
                  if (strand == PLUS) {
                     sprintf(&alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-1\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, r->ref_begin1+1+posInContigPrimer, r->ref_end1+1+posInContigPrimer, contig_length, print_strand );
                  } else {
                     sprintf(&alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-1\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, contig_length-(r->ref_end1+posInContigPrimer), contig_length - (r->ref_begin1+posInContigPrimer), contig_length, print_strand );
                  }
               }
               
               found = 1;
               break;
            }
            if (read_alignments[ind] == -1 ) {
               break;
            }
         }
         
         if (found == 0) {
            assert(bound < maxAlignments);
            read_alignments[bound] = front_contig_id;
            read_alignments_scores[bound] = r->score1;
            print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
            
            /* Translate coordinates to parent contig */
            contig_length = local_contig.parentLength;
            if (strand == PLUS) {
               posInContigPrimer = posInContigPrimer + local_contig.offsetInParent;
            } else {
               posInContigPrimer = posInContigPrimer + contig_length - (local_contig.length + local_contig.offsetInParent);
            }
            front_contig_id = local_contig.parentContigID;
            
            if (strand == PLUS) {
               sprintf(&alignment_buffer[bound * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-1\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, r->ref_begin1+1+posInContigPrimer, r->ref_end1+1+posInContigPrimer, contig_length, print_strand );
            } else {
               sprintf(&alignment_buffer[bound * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-1\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, contig_length-(r->ref_end1+posInContigPrimer), contig_length - (r->ref_begin1+posInContigPrimer), contig_length, print_strand );
            }
            bound++;
         }
         
         matched_pos_read = r->read_end1;
         align_destroy(r);
         
      }
      free(localContigInfoArray); localContigInfoArray = NULL;
   }
   
   front_ptr++;
   
   /* Try to match as much we can from this read */
   
   while ( front_ptr < (read_length - KMER_LENGTH +1) ) {
      lookup_res = NULL;
      prev_front_contig_id = front_contig_id;
      prev_front_strand = strand;

      /* Lookup for the first "valid" kmer in the REMAINING read */
      foundSeed = 0;
      while ( (lookup_res == NULL) && (front_ptr+KMER_LENGTH -1 < read_length ) ) {
         memcpy(cur_kmer, cur_read+front_ptr, KMER_LENGTH * sizeof(char));
         if (strchr(cur_kmer,'N') != NULL) {
            front_ptr++;
         } else {
			reverseComplementKmer(cur_kmer, rc_cur_kmer);
            lex_ind = strcmp(cur_kmer, rc_cur_kmer);
            kmer_to_search = ( lex_ind < 0) ? cur_kmer : rc_cur_kmer;
#ifdef PROFILE
            start_timer = UPC_TICKS_NOW();
#endif
            
#ifdef USE_SWCACHE
            lookup_res = lookupKmerInCache(hashtable, (unsigned char*) kmer_to_search, cacheCapacity, cacheTable, cacheFlags);
#else
            lookup_res = lookup_kmer(hashtable, (unsigned char*) kmer_to_search);
#endif
            
#ifdef PROFILE
            end_timer = UPC_TICKS_NOW();
            hash_table_lookups += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
            if (lookup_res == NULL) {
               front_ptr++;
            } else {
#ifdef PROFILE
               start_timer = UPC_TICKS_NOW();
#endif
               seedEntry = (list_t) *lookup_res;
               foundSeed = 1;
#ifdef PROFILE
               end_timer = UPC_TICKS_NOW();
               hash_table_lookups += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
            }
         }
      }
      
      if (foundSeed == 0) {
         //printf("The remaining read contains no valid kmer - forward search\n");
         /* Print result to output file before exiting */
#ifdef PROFILE
         start_timer = UPC_TICKS_NOW();
#endif
         for (ind = 0; ind < bound; ind++) {
            fprintf(resultFd, "%s", &alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE]);
         }
#ifdef PROFILE
         end_timer = UPC_TICKS_NOW();
         output_IO += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif
         free(alignment_buffer); alignment_buffer = NULL;
         free(read_alignments); read_alignments = NULL;
         free(read_alignments_scores); read_alignments_scores = NULL;
         assert(num);
         free(num); num = NULL;
         assert(profile);
         init_destroy(profile); profile = NULL;
         
         return bound;
      }
      
      cur_contig = seedEntry.firstContig.my_contig;
      upc_memget(&local_contig, cur_contig, sizeof(contig_t));
      contig_length = local_contig.length;
      front_contig_id = local_contig.contig_id;
      contig_copy = resetRawBuffer(contig_copy_buf, contig_length);
      assert(contig_copy[contig_length] == '\0');
      posInContig = seedEntry.firstContig.posInContig;
      assert(posInContig < contig_length);
#ifdef PROFILE
      start_fetching_timer = UPC_TICKS_NOW();
#endif
      
#ifdef USE_SWCACHE
      remote_contig_ptr = findRemoteContigPtr(front_contig_id, contig_length, cur_contig, 0, cachePtr, cacheTableContig);
#else
      remote_contig_ptr = (contigDataPtr) (&(cur_contig->contig[0]));
#endif
#ifdef DEBUG2
      fprintf(logFD, "Thread %d: Copying %d seq for contig %d %ld\n", MYTHREAD, contig_length, front_contig_id, local_contig.parentContigID);
#endif
      upc_memget(contig_copy, remote_contig_ptr, contig_length * sizeof(char));

#ifdef PROFILE
      end_fetching_timer = UPC_TICKS_NOW();
      fetch_contigs_time += (UPC_TICKS_TO_SECS(end_fetching_timer-start_fetching_timer));
#endif
      assert( contig_copy[contig_length] == '\0' );
      assert(posInContig >= 0 && contig_length >= posInContig + KMER_LENGTH);
      
      /* Start searching the position that the read and the contig start matching */
      strand = PLUS;
      if ( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) != 0 ) {
         strand =  MINUS;
#ifdef DEBUG
         //fprintf(stderr, "Thread %d: revcomp4 contig %d (%d len): %.*s\n", MYTHREAD, front_contig_id, contig_length, contig_length, contig_copy);
#endif
         reverseComplementINPLACE(contig_copy, contig_length);
         posInContig = contig_length - (KMER_LENGTH + posInContig);
         assert(posInContig < contig_length);
         assert(posInContig >= 0);
         assert( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) == 0);
      }
      
      front_contig_id = local_contig.parentContigID;
      if ((prev_front_contig_id == front_contig_id) && (prev_front_strand == strand)) {
         contig_copy = NULL;
      } else {
         /* Use SSW to align as much as we can from this read */
         /* Allow slack equal to read length since it is not that bad... */
         ref_num = (int8_t*) malloc_chk(contig_length);
         posInContigPrimer = ( (posInContig - front_ptr - read_length) < 0 ) ? 0 : posInContig - front_ptr - read_length;
         assert(posInContigPrimer < contig_length);
         assert(posInContigPrimer >= 0);
         ref_seq = contig_copy + posInContigPrimer;
         effective_length = contig_length - posInContigPrimer;
         for (m = 0; m < effective_length; ++m) ref_num[m] = (int8_t) nt_table[(int)ref_seq[m]];
         min_len = (read_length/2 >= 15) ? read_length / 2 : 15;
         assert(profile);
         assert(ref_num);
         r = ssw_align (profile, ref_num, effective_length, gap_open, gap_extension, 8, 0, 0, min_len);
         free(ref_num); ref_num = NULL;
         found = 0;
         for (ind = 0; ind < bound; ind++) {
            assert(ind < maxAlignments);
            /* Search if alignment with the same contig has been already found and check scores */
            if (read_alignments[ind] == front_contig_id) {
               /* Already have found an alignment with the same contig, check if the score is higer now and if so update output result */
               if (r->score1 > read_alignments_scores[ind]) {
                  read_alignments_scores[ind] = r->score1;
                  print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
                  
                  /* Translate coordinates to parent contig */
                  contig_length = local_contig.parentLength;
                  if (strand == PLUS) {
                     posInContigPrimer = posInContigPrimer + local_contig.offsetInParent;
                  } else {
                     posInContigPrimer = posInContigPrimer + contig_length - (local_contig.length + local_contig.offsetInParent);
                  }
                  front_contig_id = local_contig.parentContigID;
                  
                  if (strand == PLUS) {
                     sprintf(&alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, r->ref_begin1+1+posInContigPrimer, r->ref_end1+1+posInContigPrimer, contig_length, print_strand );
                  } else {
                     sprintf(&alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, contig_length-(r->ref_end1+posInContigPrimer), contig_length - (r->ref_begin1+posInContigPrimer), contig_length, print_strand );
                  }
               }
            
               found = 1;
               break;
            }
            if (read_alignments[ind] == -1 ) {
               break;
            }
         }
      
         if (found == 0) {
            assert(bound < maxAlignments); 
            read_alignments[bound] = front_contig_id;
            read_alignments_scores[bound] = r->score1;
            print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
            
            /* Translate coordinates to parent contig */
            contig_length = local_contig.parentLength;
            if (strand == PLUS) {
               posInContigPrimer = posInContigPrimer + local_contig.offsetInParent;
            } else {
               posInContigPrimer = posInContigPrimer + contig_length - (local_contig.length + local_contig.offsetInParent);
            }
            front_contig_id = local_contig.parentContigID;
            
            if (strand == PLUS) {
               sprintf(&alignment_buffer[bound * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, r->ref_begin1+1+posInContigPrimer, r->ref_end1+1+posInContigPrimer, contig_length, print_strand );
            } else {
               sprintf(&alignment_buffer[bound * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, contig_length-(r->ref_end1+posInContigPrimer), contig_length - (r->ref_begin1+posInContigPrimer), contig_length, print_strand );
            }
            bound++;
         }
      
         align_destroy(r);
      }
      
      /* If we have multiple contigs that this seed is extracted from, find alignments with those as well */
      if ( seedEntry.nExtraContigs > 0 ) {
         if (localContigInfoArray) free(localContigInfoArray);
         localContigInfoArray = (contigInfo*) malloc_chk(seedEntry.nExtraContigs * sizeof(contigInfo));
         assert( seedEntry.contigInfoArrayPtr != NULL );
#ifdef DEBUG
         fprintf(logFD, "Thread %d: getting(2) %lld contigInfoArrayPointers for %s from %d\n", MYTHREAD, (long long) seedEntry.nExtraContigs, cur_kmer, (int) upc_threadof(seedEntry.contigInfoArrayPtr));
#endif
         assert( IS_VALID_UPC_PTR( seedEntry.contigInfoArrayPtr ) );
         upc_memget(localContigInfoArray, seedEntry.contigInfoArrayPtr, seedEntry.nExtraContigs * sizeof(contigInfo));
         for (z = 0; z < seedEntry.nExtraContigs; z++) {
            /* Copy remote contig to manipulate locally */
            cur_contig = localContigInfoArray[z].my_contig;
            upc_memget(&local_contig, cur_contig, sizeof(contig_t));
            contig_length = local_contig.length;
            contig_copy = resetRawBuffer(contig_copy_buf, contig_length);
            assert(contig_copy[contig_length] == '\0');
            front_contig_id = local_contig.contig_id;
            
#ifdef PROFILE
            start_fetching_timer = UPC_TICKS_NOW();
#endif
            
#ifdef USE_SWCACHE
            remote_contig_ptr = findRemoteContigPtr(front_contig_id, contig_length, cur_contig, 0, cachePtr, cacheTableContig);
#else
            remote_contig_ptr = (contigDataPtr) (&(cur_contig->contig[0]));
#endif
#ifdef DEBUG2
            fprintf(logFD, "Thread %d: Copying %d seq for contig %d %ld\n", MYTHREAD, contig_length, front_contig_id, local_contig.parentContigID);
#endif
            upc_memget(contig_copy, remote_contig_ptr, contig_length * sizeof(char));
            
#ifdef PROFILE
            end_fetching_timer = UPC_TICKS_NOW();
            fetch_contigs_time += (UPC_TICKS_TO_SECS(end_fetching_timer-start_fetching_timer));
#endif
            posInContig = localContigInfoArray[z].posInContig;
            front_contig_id = local_contig.parentContigID;
            assert(posInContig >= 0 && contig_length >= posInContig + KMER_LENGTH); 
            
            /* Check if we should align with respect to the current contig or the reverse complement */
            strand = PLUS;
            if ( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) != 0 ) {
               strand =  MINUS;
#ifdef DEBUG
               //fprintf(stderr, "Thread %d: revcomp contig %d (%d len): %.*s\n", MYTHREAD, front_contig_id, contig_length, contig_length, contig_copy);
#endif
               reverseComplementINPLACE(contig_copy, contig_length);
               posInContig = contig_length - (KMER_LENGTH + posInContig);
               assert( memcmp(cur_kmer, contig_copy+posInContig, KMER_LENGTH) == 0);
            }
            
            /* Use SSW to align as much as we can from this read (align from start of the read to allow misalignments at the beginning)*/
            assert(!ref_num);
            ref_num = (int8_t*) malloc_chk(contig_length * sizeof(int8_t));
            read_seq = cur_read;
            posInContigPrimer = ( (posInContig - front_ptr - read_length) < 0 ) ? 0 : posInContig - front_ptr - read_length;
            assert(posInContigPrimer < contig_length);
            ref_seq = contig_copy + posInContigPrimer;
            effective_length = contig_length - posInContigPrimer;
            for (m = 0; m < effective_length; ++m) ref_num[m] = (int8_t) nt_table[(int)ref_seq[m]];
            min_len = (read_length/2 >= 15) ? read_length / 2 : 15;
            assert(profile);
            assert(ref_num);
            r = ssw_align (profile, ref_num, effective_length, gap_open, gap_extension, 8, 0, 0, min_len);
            
            free(ref_num); ref_num = NULL;
            found = 0;
            for (ind = 0; ind < bound; ind++) {
               assert(ind < maxAlignments);
               /* Search if alignment with the same contig has been already found and check scores */
               if (read_alignments[ind] == front_contig_id) {
                  /* Already have found an alignment with the same contig, check if the score is higer now and if so update output result */
                  if (r->score1 > read_alignments_scores[ind]) {
                     read_alignments_scores[ind] = r->score1;
                     print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
                     
                     /* Translate coordinates to parent contig */
                     contig_length = local_contig.parentLength;
                     if (strand == PLUS) {
                        posInContigPrimer = posInContigPrimer + local_contig.offsetInParent;
                     } else {
                        posInContigPrimer = posInContigPrimer + contig_length - (local_contig.length + local_contig.offsetInParent);
                     }
                     front_contig_id = local_contig.parentContigID;
                     
                     if (strand == PLUS) {
                        sprintf(&alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, r->ref_begin1+1+posInContigPrimer, r->ref_end1+1+posInContigPrimer, contig_length, print_strand );
                     } else {
                        sprintf(&alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, contig_length-(r->ref_end1+posInContigPrimer), contig_length - (r->ref_begin1+posInContigPrimer), contig_length, print_strand );
                     }
                  }
                  
                  found = 1;
                  break;
               }
               if (read_alignments[ind] == -1 ) {
                  break;
               }
            }
            
            if (found == 0) {
               assert(bound < maxAlignments);
               read_alignments[bound] = front_contig_id;
               read_alignments_scores[bound] = r->score1;
               print_strand = ( strand == PLUS) ? plus_strand : minus_strand;
               
               /* Translate coordinates to parent contig */
               contig_length = local_contig.parentLength;
               if (strand == PLUS) {
                  posInContigPrimer = posInContigPrimer + local_contig.offsetInParent;
               } else {
                  posInContigPrimer = posInContigPrimer + contig_length - (local_contig.length + local_contig.offsetInParent);
               }
               front_contig_id = local_contig.parentContigID;
               
               if (strand == PLUS) {
                  sprintf(&alignment_buffer[bound * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, r->ref_begin1+1+posInContigPrimer, r->ref_end1+1+posInContigPrimer, contig_length, print_strand );
               } else {
                  sprintf(&alignment_buffer[bound * MAX_ALIGNMENT_LINE_SIZE], "MERALIGNER-2\t%s\t%d\t%d\t%d\tContig%d\t%d\t%d\t%d\t%s\t0\t0\t0\t0\n", read_info, r->read_begin1+1, r->read_end1+1, read_length, front_contig_id, contig_length-(r->ref_end1+posInContigPrimer), contig_length - (r->ref_begin1+posInContigPrimer), contig_length, print_strand );
               }
               bound++;
            }
            
            matched_pos_read = r->read_end1;
            align_destroy(r);
            
         }
         assert(localContigInfoArray);
         free(localContigInfoArray); localContigInfoArray = NULL;
      }
   
      front_ptr++;
   }
   
   /* Print result to output file before exiting */
#ifdef PROFILE
   start_timer = UPC_TICKS_NOW();
#endif
   for (ind = 0; ind < bound; ind++) {
      fprintf(resultFd, "%s", &alignment_buffer[ind * MAX_ALIGNMENT_LINE_SIZE]);
   }
#ifdef PROFILE
   end_timer = UPC_TICKS_NOW();
   output_IO += (UPC_TICKS_TO_SECS(end_timer-start_timer));
#endif

   assert(alignment_buffer);
   assert(read_alignments);
   assert(read_alignments_scores); 
   free(alignment_buffer); alignment_buffer = NULL;
   free(read_alignments); read_alignments = NULL;
   free(read_alignments_scores); read_alignments_scores = NULL;

   assert(profile);
   init_destroy(profile); profile = NULL;
   assert(!ref_num);
   if (num) free(num); num = NULL;

   return bound;
}


#endif

