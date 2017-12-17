#ifndef __COMMON_H
#define __COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <sys/stat.h>
#include <stdint.h>
#include <ctype.h>
#include <assert.h>
#include <stdarg.h>
#include <unistd.h>

#include "version.h"

#include "colors.h"

#define MAX_FILE_PATH 384

#ifdef USE_UPC_FOR_COMMON

#include <upc.h>
#define EXIT_FUNC(x) do { upc_global_exit(x); } while(0)

#elif defined USE_MPI_FOR_COMMON

#include <mpi.h>
#define EXIT_FUNC(x) do { MPI_Abort(MPI_COMM_WORLD, x); exit(x); } while (0)
static int get_rank(void) {int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank); return rank;}
static int get_num_ranks(void) {int num_ranks; MPI_Comm_size(MPI_COMM_WORLD, &num_ranks); return num_ranks;}
#define MYTHREAD get_rank()
#define THREADS get_num_ranks()

#else

#error Must define either USE_UPC_FOR_COMMON or USE_MPI_FOR_COMMON

#endif 

#define DIE(fmt,...)													\
    do {                                                                \
        fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: " fmt KNORM,		\
                MYTHREAD, __FILE__, __LINE__, HIPMER_VERSION, ##__VA_ARGS__);           \
        fflush(stderr);                                                \
        EXIT_FUNC(1);                                                   \
    } while (0)

#define SDIE(fmt,...)													\
    do {                                                                \
        if (!MYTHREAD)                                                  \
            fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: " fmt KNORM, \
                    MYTHREAD, __FILE__, __LINE__, HIPMER_VERSION, ##__VA_ARGS__);       \
        fflush(stderr);                                                \
        EXIT_FUNC(1);                                                   \
    } while (0)

#define WARN(fmt,...)                                               \
    do {                                                            \
        fprintf(stderr, KRED "Thread %d, WARN [%s:%d-%s]: " fmt KNORM, \
                MYTHREAD, __FILE__, __LINE__, HIPMER_VERSION, ##__VA_ARGS__);       \
    } while (0)

#define serial_printf(...)                      \
    do {                                        \
        if (MYTHREAD == 0) {                    \
            fprintf(stderr, __VA_ARGS__);       \
            fflush(stderr);                     \
        }                                       \
    } while (0)

#define CHECK_ERR(cmd)                                                  \
    do {                                                                \
        int err;                                                        \
        if ((err = cmd) != 0)                                           \
            DIE("Thread %d, [%s:%d-%s]: " #cmd " failed, error %s\n", MYTHREAD, __FILE__, __LINE__, HIPMER_VERSION, strerror(err)); \
    } while (0)

#define INT_CEIL(numerator, denominator) (((numerator) - 1) / (denominator) + 1);

#define sprintf_chk(s, max, fmt, ...) \
    do { \
        int linelen = sprintf(s, fmt, ##__VA_ARGS__); \
        if (linelen > max) \
           DIE("Thread %d: buffer overflow printing '%s' at %s:%d-%s: (" fmt ")\n", MYTHREAD, fmt,  __FILE__, __LINE__, HIPMER_VERSION, ##__VA_ARGS__); \
    } while (0)


inline off_t get_file_size(const char *fname)
{
    struct stat s;
    if (stat(fname, &s) != 0)
        DIE("could not stat %s: %s\n", fname, strerror(errno));
    return s.st_size;
}

#define malloc_chk(s) malloc_chk_at_line(s, __FILE__, __LINE__)
inline void *malloc_chk_at_line(size_t size, const char *which_file, int which_line)
{
	if (size == 0) return NULL;
	void *x = malloc(size);
	if (!x) {
		fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: Could not malloc %lu bytes\n" KNORM, 
			MYTHREAD, which_file, which_line, HIPMER_VERSION, size);
		EXIT_FUNC(1);
	}
	return x;
}

#define calloc_chk(n, s) calloc_chk_at_line(n, s, __FILE__, __LINE__)
inline void *calloc_chk_at_line(size_t nmemb, size_t size, const char *which_file, int which_line)
{
	if (size == 0) return NULL;
	void *x = calloc(nmemb, size);
	if (!x) {
		fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: Could not calloc %lu bytes\n" KNORM, 
			MYTHREAD, which_file, which_line, HIPMER_VERSION, size);
		EXIT_FUNC(1);
	}
	return x;
}

#define realloc_chk(p, s) realloc_chk_at_line(p, s, __FILE__, __LINE__)
inline void *realloc_chk_at_line(void *ptr, size_t size, const char *which_file, int which_line)
{
	assert(size>0);
	void *x = realloc(ptr, size);
	if (!x) {
		fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: Could not realloc %lu bytes\n" KNORM, 
			MYTHREAD, which_file, which_line, HIPMER_VERSION, size);
		EXIT_FUNC(1);
	}
	return x;
}

#define fopen_chk(p, m) fopen_chk_at_line(p, m, __FILE__, __LINE__)
inline FILE *fopen_chk_at_line(const char *path, const char *mode, const char *which_file, int which_line)
{
    FILE *f = fopen(path, mode);
    if (!f) {
        fprintf(stderr, KLRED "Thread %d, DIE [%s:%d-%s]: Could not open file '%s': %s\n" KNORM, 
                MYTHREAD, which_file, which_line, HIPMER_VERSION, path, strerror(errno));
        EXIT_FUNC(1);
    }
    return f;
}

inline char *get_basename(char *bname, const char *path)
{
	const char *slash_pos = rindex(path, '/');
	if (!slash_pos) 
		strcpy(bname, path);
	else 
		strcpy(bname, slash_pos+1);
	return bname;
}

inline char *get_dirname(char *dname, const char *path)
{
	const char *slash_pos = rindex(path, '/');
	if (!slash_pos) 
		strcpy(dname, path);
	else 
		strncpy(dname, path, (slash_pos - path));
	return dname;
}

inline int check_seqs(const char *seq, char *label)
{
    int len = strlen(seq);
	for (int i = 0; i < len; i++) {
        char c = toupper(seq[i]);
		if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') 
			DIE("Invalid base %c %d in sequence (len %d) %.80s\n%s\n", c, c, len, seq, label);
	}
    return len;
}

// returns 1 when it created the directory, 0 otherwise, dies if there is an error
inline int check_dir(const char *path) {
    if (0 != access(path, F_OK)) {
       if (ENOENT == errno) {
          // does not exist
          if (0 != mkdir(path, 0777)) {
              DIE("Could not create the (missing) directory: %s (%s)", path, strerror(errno));
          } 
       }
       if (ENOTDIR == errno) {
           // not a directory
           DIE("Expected %s was a directory!", path);
       }
    } else {
       return 0;
    }
    assert( access(path, F_OK) == 0 );
    return 1;
}

// replaces the given path with a rank based path, inserting a rank-based directory
// example:  get_rank_path("/path/to/file_output_data.txt", rank) -> "path/to/per_thread/<rank>/file_output_data.txt"
// of if rank == -1, "path/to/per_thread/file_output_data.txt"
inline char * get_rank_path(char *buf, int rank) {
    int pathlen = strlen(buf);
    if (pathlen+15 >= MAX_FILE_PATH) {
       DIE("File path is too long (max: %d): %s\n", MAX_FILE_PATH, buf);
    }

    // get the last basedir/basename
    int padding = rank < 0 ? 0 : 8;
    const char *prefix = "per_thread/";
    int prefixlen = strlen(prefix);
    char *slash_pos = rindex(buf, '/');
    int slash_idx = 0;
    int move_start = 0;
    if (slash_pos) {
      // there was a basedir: "basedir/file"
      assert(*slash_pos == '/');
      slash_idx = (int) (slash_pos - buf);
      assert(buf[slash_idx] == '/');
      move_start = slash_idx + 1; // move the directory, leave the slash
    }
    memmove(buf + move_start + padding + prefixlen, buf + move_start, pathlen - move_start);
    memcpy(buf+move_start, prefix, prefixlen);
    move_start += prefixlen;
    pathlen += prefixlen;
    slash_pos = buf + move_start -1;
    slash_idx = (int) (slash_pos - buf);
    assert(buf[slash_idx] == '/');
    buf[slash_idx] = '\0';
    check_dir(buf);
    buf[slash_idx] = '/';
    assert(buf[slash_idx] == '/');
 
    // nopath: slash_idx=0,move_start=0,padding=8: "        file"
    // path:   slash_idx=7,move_start=8,padding=8: "basedir/        file"
    
    // ensure all strings and substrings are terminated
    buf[pathlen + padding] = '\0';
    // "|       file|"
    // "basedir/|       file|"

    // insert and create the rankdir
    if (rank >= 0) {
      int tmp = move_start + padding -1;
      buf[tmp] = '\0';
      int printrank = rank;
      while (tmp > move_start) {
          buf[--tmp] = (char) ('0' + (printrank%10));
          printrank /= 10;
      }
      // "7654321|file|"
      // "basedir/7654321|file|"
    
      assert(tmp == move_start);
      assert(buf[tmp] == '0');
      assert(buf[slash_idx] == '/');
      assert(strlen(buf) == move_start + padding - 1);
      int wasCreated = check_dir(buf);

      if (wasCreated) {
          // pre-emptively set lustre striping to 1 (do not fail if this filesystem is not lustre
          char cmd[MAX_FILE_PATH+35];
          sprintf(cmd, "lfs setstripe -c 1 %s 2>/dev/null", buf);
          if (system(cmd) == -1) 
              WARN("Error executing %s\n", cmd);
      }
      assert(buf[move_start + padding - 1] == '\0');
    } else {
      assert(padding==0);
      assert(buf[move_start + padding - 1] == '/');
    }
    buf[move_start + padding - 1] = '/';
    // "7654321/file|"
    // "basedir/7654321/file|"

#ifdef DEBUG
    fprintf(stderr, "Thread %d: converted to rank_path '%s'\n", rank, buf);
#endif
    return buf;
}

#define fopen_rank_path(p, m, r) fopen_rank_path_at_line(p, m, r, __FILE__, __LINE__)
static FILE * fopen_rank_path_at_line(char *buf, const char *mode, int rank, const char *which_file, int which_line) {
    char *rankpath = get_rank_path(buf, rank);
    FILE *f = fopen_chk_at_line(rankpath, mode, which_file, which_line);
    return f;
}

#if defined (__cplusplus)
#include <string>
static std::string getRankPath(const char *path, int rank) {
    char buf[MAX_FILE_PATH];
    strcpy(buf, path);
    get_rank_path(buf, rank);
    return std::string(buf);
}
#endif

#define UPC_MEMGET_STR(dst, src, len) do {                              \
        upc_memget(dst, src, len);                                      \
        if (strlen(dst) != len - 1)                                     \
            DIE("Length mismatch when getting remote NULL terminated string: %ld != %d\n", strlen(dst), len - 1); \
    } while (0)

#ifndef SEGMENT_LENGTH
#define SEGMENT_LENGTH 51
#endif

static size_t printFoldedSequence(FILE *out, const char *finalSequence, size_t cur_length) {
    /* Print contig in FASTA FORMAT */
    size_t total_written = 0, towrite;
    const char *seqF = (char*) finalSequence;
    char fastaSegment[SEGMENT_LENGTH];
    while ( total_written < cur_length ) {
        if (total_written + SEGMENT_LENGTH-1 < cur_length) {
            towrite = SEGMENT_LENGTH-1;
        } else {
            towrite = cur_length - total_written;
        }
        memcpy(fastaSegment, seqF+total_written, towrite * sizeof(char));
        fastaSegment[towrite] = '\0';
        fprintf(out, "%s\n", fastaSegment);
        total_written += towrite;
    }
    assert(total_written == cur_length);
    return total_written;
}

uint32_t MurmurHash3_x64_32(const void *key, int len);
uint64_t MurmurHash3_x64_64(const void *key, int len);
uint64_t murmur_hash2_64(const void * key, int len);

//#define hashkey(table_size, key, len) (MurmurHash3_x64_64(key, len) % (table_size))
//#define hashstr(table_size, key) (MurmurHash3_x64_64(key, strlen(key)) % (table_size))
#define hashkey(table_size, key, len) (murmur_hash2_64(key, len) % (table_size))
#define hashstr(table_size, key) (murmur_hash2_64(key, strlen(key)) % (table_size))

#endif
