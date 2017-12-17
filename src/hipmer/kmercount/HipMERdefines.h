#ifndef __HIPMERDEFINES_H
#define __HIPMERDEFINES_H

/*
#if defined ILLUMINA_VERSION && ILLUMINA_VERSION >= 10 && ILLUMINA_VERSION < 15
        #define PHREDENCODING 64        // Phred+64
#elif defined ILLUMINA_VERSION && ILLUMINA_VERSION >= 15 && ILLUMINA_VERSION < 18
        #define PHREDENCODING 64        // Phred+64
#elif defined ILLUMINA_VERSION && ILLUMINA_VERSION >= 18
        #define PHREDENCODING 33        // Phred+33
#else
        #define PHREDENCODING 33        // default
#endif
*/
// Always trim the read at qual==2
#define TRIMSPECIAL // http://en.wikipedia.org/wiki/FASTQ_format#Encoding       

#ifndef KMER_LENGTH
        #define KMER_LENGTH 51
#endif

#ifndef MAX_KMER_SIZE
        #if defined KMER_LENGTH && KMER_LENGTH <= 32
                #define MAX_KMER_SIZE 32
        #elif defined KMER_LENGTH && KMER_LENGTH <= 64
                #define MAX_KMER_SIZE 64
        #elif defined KMER_LENGTH && KMER_LENGTH <= 96
                #define MAX_KMER_SIZE 96
        #elif defined KMER_LENGTH && KMER_LENGTH <= 128
                #define MAX_KMER_SIZE 128
        #else
                #define MAX_KMER_SIZE 160
        #endif
#endif

#define KMERLONGS MAX_KMER_SIZE/32      // 32 = numbits(uint64_t)/2-  with 2 being the number of bits needed per nucleotide

static int phred_encoding = 33;

#define extensionCutoff 20

static int ignored_ext = phred_encoding + 1;// too low to be a good extension
//#define IGNOREDEXT PHREDENCODING+1 

#endif
