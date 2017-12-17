/*
 * Note: this fastq reader code is not specific to upc. It is also used by kMerCount with mpi.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>

#include "fq_reader.h"

#define NUM_STRIDES 100

#ifndef MMAP_BLOCK_SIZE
#define MMAP_BLOCK_SIZE getpagesize()
#endif

fq_reader_t create_fq_reader(void)
{ 
    fq_reader_t fqr = (fq_reader_t) calloc_chk(1, sizeof(struct fq_reader));
    fqr->buf = initBuffer(MAX_READ_LEN);
    assert(isValidBuffer(fqr->buf));
    return fqr;
}

void destroy_fq_reader(fq_reader_t fqr)
{
    if (fqr) {
        if (fqr->buf) 
            freeBuffer(fqr->buf);
        fqr->buf = NULL;

        if (fqr->f) 
            close_fq(fqr);
        fqr->f = NULL;

        if (fqr->addr && fqr->size > 0) { 
            size_t offset = fqr->start_read % MMAP_BLOCK_SIZE;
            munmap(fqr->addr, fqr->end_read - fqr->start_read + offset);
        }
        fqr->addr = NULL;
        fqr->size = 0;
        free(fqr);
    }
}


static void strip_trailing_spaces(char *s)
{
    char *end = s + strlen(s) - 1;
    while (end >= s && isspace(*end))
        end--;
    *(end + 1) = '\0';
}

inline static char *get_fq_name(char *header)
{
    assert(header != NULL);
    assert(header[0] == '@');
    int len = strlen(header);
    if (len >= MAX_READ_NAME_LEN - 1) {
        return NULL;
    }
    strip_trailing_spaces(header);
    len = strlen(header);
    // only convert if new illumina 2 format
    if (header[len - 2] != '/') {
        // Latest Illumina header format
        char *end = strchr(header, '\t');
        if (!end) {
            end = strchr(header, ' ');
            if (!end) {
                // no comment, just return the name without modification
                return header;
            }
        }
        // truncate name at comment
        assert(*end == ' ' || *end == '\t');
        end[0] = '\0';
        // check for @pair/1 @/pair2 vs @pair 1:Y:... @pair 2:Y:....
        int commentpos = end - header;
        if (commentpos > 3 && header[commentpos-2] == '/' && (header[commentpos-1] == '1' || header[commentpos-1] == '2')) {
            // @pair/1 or @pair/2
            return header;
        } else if ((len < commentpos + 7) || end[2] != ':' || end[4] != ':' || end[6] != ':' || (end[1] != '1' && end[1] != '2')) {
            // unknown pairing format
            return header;
        }
        // @pair 1:Y:.:.: or @pair 2:Y:.:.:
        // replace with @pair/1 or @pair/2
        end[0] = '/';
        end[2] = '\0';
    }
    return header;
}

int get_fq_name_dirn(char *header, char **name, int *end)
{
    *name = get_fq_name(header);
    if (!*name) {
        return 0;
    }
    int len = strlen(*name);
    if (len < 3) {
        return 0;
    }
    char end_char = (*name)[len - 1];
    if ((*name)[len-2] != '/' || (end_char != '1' && end_char != '2')) {
        return 0;
    }
    *end = (end_char == '1' ? 1 : 2);
    (*name)[len - 2] = '\0';
    return 1;
}

// returns 0 for unpaired (offset to '@'), or the offset to '1' or '2' for paired
static int get_pairidx_from_name_line(char *s)
{
    int pairIdx = 0, max = 2047;
    // Supports two types of paired formats:
    //
    // @name/1
    // @name/2  
    // 
    // @HISEQ15:136:C6G2YANXX:4:1101:1074:2037 1:N:0:AGTTCC
    // @HISEQ15:136:C6G2YANXX:4:1101:1074:2037 2:N:0:AGTTCC
    
    assert(s[0] == '@');
    int endIdx = 1, commentIdx = 0;
    while (s[endIdx] != '\n' && s[endIdx] != '\0') {
        if(endIdx >= max) break;
        if (commentIdx == 0 && isspace(s[endIdx])) commentIdx = endIdx;
        endIdx++;
    }
    if (s[endIdx] != '\n') return pairIdx; // line too long to determine pairing!
    if (commentIdx == 0) commentIdx = endIdx;
    if (commentIdx > 3 && s[commentIdx-2] == '/' && (s[commentIdx-1] == '1' || s[commentIdx-1] == '2')) {
        pairIdx = commentIdx-1;
    } else if (endIdx > commentIdx + 7 && (s[commentIdx+2] == ':' && s[commentIdx+4] == ':' && s[commentIdx+6] == ':' &&
               (s[commentIdx+1] == '1' || s[commentIdx+1] == '2'))) {
        pairIdx = commentIdx+1;
        
    }
    return pairIdx;
}

static int64_t get_fptr_for_next_record(fq_reader_t fqr, int64_t offset)
{
    fseek(fqr->f, offset, SEEK_SET);
    // clear from this offset to ensure we start at the beginning of a line
    if (!fgetsBuffer(fqr->buf, 2047, fqr->f)) 
        return ftell(fqr->f);
    int64_t new_offset = 0;
    int count = 0, last_pair_line = 0;
    char prev_pair = 0;
    int64_t prev_offset = 0;
    while (1) {
        new_offset = ftell(fqr->f);
        resetBuffer(fqr->buf);
        if (!fgetsBuffer(fqr->buf, 2047, fqr->f))
            break;
        count++;
        char *buf = getStartBuffer(fqr->buf);
        if (buf[0] == '@') {
            int pairIdx = get_pairidx_from_name_line(buf);
            // not paired! could be a quality line
            if (pairIdx == 0) 
                continue; 
            char pair = buf[pairIdx];
            if (!prev_pair) {
                prev_pair = pair;
                prev_offset = new_offset;
                last_pair_line = count;
            } else {
                if (last_pair_line + 1 == count) {
                    // last pair was actually a quality line
                    prev_pair = pair;
                    prev_offset = new_offset;
                    last_pair_line = count;
                } else if (last_pair_line + 4 == count) {
                    // not interleaved or interleaved and prev pair was 1 and this is 2
                    if (pair == prev_pair || (prev_pair == '1' && pair == '2') ) 
                        new_offset = prev_offset;
                    break;
                }
            }

         }
         if (count > 9) 
            DIE("Could not find a valid line in the fastq file, last line: %s\n", buf);
    }
    // for interleawed, make sure to the next record starts with a /1
    return new_offset;
}

static void set_read_block(fq_reader_t fqr)
{
    int64_t read_block = INT_CEIL(fqr->size, THREADS * NUM_STRIDES);
    int64_t offset = MYTHREAD + (fqr->stride * THREADS);
    fqr->start_read = read_block * offset;
    fqr->end_read = read_block * (offset + 1);
    if (offset > 0) 
        fqr->start_read = get_fptr_for_next_record(fqr, fqr->start_read);
    if (MYTHREAD == THREADS - 1 && fqr->stride == NUM_STRIDES - 1)
        fqr->end_read = fqr->size;
    else 
        fqr->end_read = get_fptr_for_next_record(fqr, fqr->end_read);
    //serial_printf("Stride %d, new start for thread 0 is %ld, line %ld\n", 
    //              fqr->stride, fqr->start_read, fqr->line);
}

static void open_fq_disk(fq_reader_t fqr, const char *fname, int64_t knownFileSize)
{
    assert(strlen(fname) < MAX_FILE_PATH);
    strcpy(fqr->name, fname);
    fqr->stride = 0;

    if (knownFileSize <= 0) {
        // avoid a stat (across all nodes) if the size is already known...
        fqr->size = get_file_size(fname);
    } else {
        fqr->size = knownFileSize;
        assert(fqr->size == get_file_size(fname));
    }
    fqr->f = fopen_chk(fname, "r");
    set_read_block(fqr);
    fqr->f = freopen(fname, "r", fqr->f);
    if (!fqr->f) DIE("Could not reopen %s! %s", fname, strerror(errno));
    int e = setvbuf(fqr->f, NULL, _IOFBF, FQ_READER_BUFFER_SIZE);
    if (e != 0) WARN("Could not setvbuf on %s to %d! %s", fname, FQ_READER_BUFFER_SIZE, strerror(errno));

#ifdef USE_MMAP_TO_READ_FQ
    // additionally mmap the file to passively assist in re-reading
    size_t offset = fqr->start_read % MMAP_BLOCK_SIZE, maplen = fqr->end_read - fqr->start_read;
    fqr->addr = (char*) mmap(NULL, maplen + offset, PROT_READ, MAP_SHARED, fileno(fqr->f), 
                             fqr->start_read - offset);
    if (fqr->addr == MAP_FAILED) 
        fqr->addr = NULL;
#else
    assert(fqr->addr == NULL);
#endif
    reset_my_partition(fqr);
    assert(fqr->fpos == fqr->start_read);
    assert(fqr->fpos == ftell(fqr->f));
}

static int open_my_fq_shm(fq_reader_t fqr, const char *fname, int64_t knownFileSize)
{
    char *my_fname = fqr->name;
    char bname[MAX_FILE_PATH];
    sprintf(my_fname, "%s.%d", get_basename(bname, fname), MYTHREAD);
    int fd = shm_open(my_fname, O_RDONLY, 1);
    if (fd == -1) 
        return 0;

    if (knownFileSize <= 0) {
        struct stat sbuf;
        fstat(fd, &sbuf);
        knownFileSize = sbuf.st_size;
    }

    fqr->addr = (char*)mmap(NULL, knownFileSize, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (fqr->addr == MAP_FAILED) 
        return 0;
    fqr->size = knownFileSize;
    // for keeping track of when we are done and progress
    fqr->start_read = 0;
    fqr->end_read = knownFileSize;
    fqr->fpos = 0;
    fqr->line = 0;
    return 1;
}

void open_fq(fq_reader_t fqr, const char *fname, int from_shm) {
    open_fq_size(fqr, fname, from_shm, -1);
}
void open_fq_size(fq_reader_t fqr, const char *fname, int from_shm, int64_t knownFileSize)
{
    if (knownFileSize <= 0)
        knownFileSize = get_file_size(fname);
    if (from_shm) {
        if (!open_my_fq_shm(fqr, fname, knownFileSize))
            DIE("Could not open shared memory for %s\n", fname);
    } else {
        open_fq_disk(fqr, fname, knownFileSize);
    }
#ifdef CONFIG_SHOW_PROGRESS
    fqr->tick_size = (fqr->end_read - fqr->start_read) / 10;
    fqr->tick = 1;
    serial_printf("Reading FASTQ file %s %s: ", fname, from_shm ? " (shm)" : "");
    fqr->done = 0;
#endif
    fqr->max_read_len = 0;
}

void reset_my_partition(fq_reader_t fqr)
{
    if (fqr->addr != NULL) {
        int64_t offsetpos = fqr->start_read % MMAP_BLOCK_SIZE;
        madvise(fqr->addr + fqr->start_read - offsetpos, 
                fqr->end_read - fqr->start_read + offsetpos, MADV_SEQUENTIAL);
    }
    if (fqr->f != NULL)
        fseek(fqr->f, fqr->start_read, SEEK_SET);
    fqr->fpos = fqr->start_read;
    fqr->line = 0;
#ifndef __APPLE__
    posix_fadvise(fileno(fqr->f), fqr->start_read, fqr->end_read - fqr->start_read, 
                  POSIX_FADV_SEQUENTIAL);
#endif
}

int load_fq(fq_reader_t fqr, char *fname)
{
    open_fq_disk(fqr, fname, -1);
    int err = 0;
    // create shm file
    char my_fname[MAX_FILE_PATH];
    char bname[MAX_FILE_PATH];
    sprintf(my_fname, "%s.%d", get_basename(bname, fname), MYTHREAD);
    int fd = shm_open(my_fname, O_CREAT|O_TRUNC|O_RDWR, S_IRUSR|S_IWUSR);
    if (fd != -1) {
        int64_t read_block = fqr->end_read - fqr->start_read;
        if (ftruncate(fd, read_block) == -1) {
            WARN("Failed to truncate %s: %s\n", my_fname, strerror(errno));
            err = errno;
        } else {
            void *tmpaddr = mmap(NULL, read_block, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
            if (tmpaddr == MAP_FAILED) {
                WARN("Could not mmap %s for writing: %s\n", my_fname, strerror(errno));
                err = errno;
            } else {
                // read in 10 blocks so we can get show progress
                int64_t sub_block = INT_CEIL(read_block, 10);
                size_t bytes_read = 0;
                for (int i = 0; i < 10; i++) {
					int64_t pos = i * sub_block;
					if (bytes_read + sub_block > read_block)
						sub_block = read_block - bytes_read;
                    bytes_read += fread((char*)tmpaddr + pos, 1, sub_block, fqr->f);
#ifdef CONFIG_SHOW_PROGRESS
                    serial_printf("%d ", i + 1);
#endif
                }
#ifdef CONFIG_SHOW_PROGRESS
                serial_printf("\n");
#endif
				//fprintf(stderr, "thread %d done with reading\n", MYTHREAD);
                if (msync(tmpaddr, bytes_read, MS_SYNC) == -1)
					WARN("msync returned error %s\n", strerror(errno));
                if (munmap(tmpaddr, bytes_read) == -1)
					WARN("munmap returned error %s\n", strerror(errno));
                if (fseek(fqr->f, fqr->start_read, SEEK_SET) == -1)
					WARN("fseek returned error %s\n", strerror(errno));
            }
        }
        close(fd);
    }
    close_fq(fqr);
    return -err;
}

int unload_fq(char *fname)
{
    char my_fname[MAX_FILE_PATH];
    char bname[MAX_FILE_PATH];
    sprintf(my_fname, "%s.%d", get_basename(bname, fname), MYTHREAD);
    if (shm_unlink(my_fname) == -1) 
        return -errno;
    return 0;
}

// include newline char
static char *get_next_line(fq_reader_t fqr)
{
    char *buf = NULL;
    size_t len = 0, maxLen = (fqr->end_read - fqr->fpos);
    resetBuffer(fqr->buf);
    if (fqr->addr) {
        char *start = fqr->addr + fqr->fpos;
        char *end = (char*) memchr(start, '\n', maxLen);
        if (!end) 
            return NULL;
        assert(*end == '\n');
        len = (end - start) + 1; // include newline
        assert(len <= maxLen);
        growBuffer(fqr->buf, len+1);
        char *buf = (char*) memcpyBuffer(fqr->buf, start, len);

        assert(getSizeBuffer(fqr->buf) > len); // has room for zero terminated string
        buf[len] = '\0';
        fqr->fpos += len;
    } else {
        if (fqr->fpos >= fqr->end_read) return NULL;
        buf = fgetsBuffer(fqr->buf, (maxLen > 2048 ? 2048 : maxLen), fqr->f);
        if (!buf)
           return NULL;
        
        len = getLengthBuffer(fqr->buf);
        fqr->fpos += len;
        assert(fqr->fpos == ftell(fqr->f));
    }
    assert(buf[len] == '\0');
    if (buf[len-1] != '\n')
        WARN("should be end of line, not '%c' %d: '%s' at line %ld\n", buf[len], buf[len], buf, fqr->line);
    fqr->line++;
    assert(buf == getStartBuffer(fqr->buf));
    assert(len > 0);
    assert(len == getLengthBuffer(fqr->buf));
    return buf;
}

int get_next_fq_record(fq_reader_t fqr, Buffer id, Buffer nts, Buffer quals) {
    resetBuffer(id);
    resetBuffer(nts);
    resetBuffer(quals);
#ifdef CONFIG_SHOW_PROGRESS
    // file progress
    if ((fqr->fpos - fqr->start_read) / fqr->tick_size == fqr->tick) {
        serial_printf("%d ", fqr->tick);
        fqr->tick++;
    }
#endif
    if (fqr->fpos >= fqr->end_read) {
#ifdef CONFIG_SHOW_PROGRESS
        if (!fqr->done) {
            serial_printf(" FASTQ done\n");
            fqr->done = 1;
        }
#endif
        fqr->stride++;
        if (fqr->stride == NUM_STRIDES)
            return 0;
        set_read_block(fqr);
        reset_my_partition(fqr);
        assert(fqr->fpos == fqr->start_read);
        assert(fqr->fpos == ftell(fqr->f));
    }
    Buffer lastBuf = fqr->buf;

    // get all four lines, one for each field
    for (int i = 0; i < 4; i++) {
        char *tmp = NULL;
        switch(i) {
            case 0: fqr->buf = lastBuf; break;
            case 1: fqr->buf = nts;     break; 
            case 2: fqr->buf = lastBuf; break;
            case 3: fqr->buf = quals;   break;
            default: assert(0);
        };
        tmp = get_next_line(fqr);
        if (!tmp) {
#ifdef CONFIG_SHOW_PROGRESS
            if (!fqr->done) {
                serial_printf(" FASTQ done\n");
                fqr->done = 1;
            }
#endif
            // return the proper buffer back
            fqr->buf = lastBuf;
            return 0;
        }

        char *buf = getStartBuffer(fqr->buf);
        if (getLengthBuffer(fqr->buf) == 0 || buf[0] == '\n') {
            // empty
            if (i != 0 ) 
                DIE("Invalid FASTQ at line %ld, expected empty or just newline at the end: '%s'\n", fqr->line, buf);
            // return the proper buffer back
            fqr->buf = lastBuf;
            return 0;
        }
            
        if (i == 0) {
            // header line
            if (buf[0] != '@')        
                DIE("Invalid FASTQ at line %ld, expected read name (@):\n%s\n", fqr->line, buf);
            // construct universally formatted name (illumina 1 format)
            strip_trailing_spaces(buf);
            char *read_name = get_fq_name(buf);
            if (!read_name)
                DIE("Invalid FASTQ name format: %s\n", buf);
            strcpyBuffer(id, read_name);
        } else if (i == 1) {
            // sequence
            chompBuffer(fqr->buf);
#ifdef CONFIG_CHECK_SEQS        
            char label[256];
            sprintf(label, "  in %s:%d at line %ld\n", __FILE__, __LINE__, fqr->line);
            assert(strlen(buf) == getLengthBuffer(fqr->buf));
            int len = check_seqs(buf, label);
#endif
            assert(fqr->buf == nts);
        } else if (i == 2) {
            if (buf [0] != '+')
                DIE("Invalid FASTQ at line %ld, expected '+':\n%s\n", fqr->line, buf);
        } else if (i == 3) {
            chompBuffer(fqr->buf);
            assert(fqr->buf == quals);
        } else {
            DIE("");
        }

    }
    if (getLengthBuffer(nts) != getLengthBuffer(quals)) {
        DIE("Sequence and quals differ in length at line %ld: %lu != %lu\n%s\n%s\n", 
            fqr->line, getLengthBuffer(nts), getLengthBuffer(quals), getStartBuffer(nts), getStartBuffer(quals));
    }
    int read_len = getLengthBuffer(nts);
    if (read_len > fqr->max_read_len)
        fqr->max_read_len = read_len;
    // reset to the original Buffer
    fqr->buf = lastBuf;
    return 1;
}

/*
int get_next_fq_record_ptr(fq_reader_t fqr, char **id, char **nts, char **quals) {
    Buffer _id = initBuffer(MAX_READ_NAME_LEN);
    Buffer _nts = initBuffer(MAX_READ_LEN);
    Buffer _quals = initBuffer(MAX_READ_LEN);
    if ( get_next_fq_record(fqr, _id, _nts, _quals) == 1) {
        if (*id) free (*id);
        if (*nts) free(*nts);
        if (*quals) free(*quals);
        *id = releaseBuffer(_id);
        *nts = releaseBuffer(_nts);
        *quals = releaseBuffer(_quals);
        return 1;
    } else {
        freeBuffer(_id);
        freeBuffer(_nts);
        freeBuffer(_quals);
        return 0;
    }
}
*/  

void close_fq(fq_reader_t fqr)
{
    if (fqr->f) {
        fclose(fqr->f);
        fqr->name[0] = '\0';
        fqr->f = NULL;
        fqr->stride = 0;
    }
}


int64_t estimate_fq(char *fname, int sampledReads, int64_t *estimatedReads, int64_t *estimatedBases)
{
    fq_reader_t fqr = create_fq_reader();
    open_fq_disk(fqr, fname, -1);
    int max = sampledReads;
    int i;
    int64_t bases = 0;
    int64_t startPos = fqr->fpos;
    for(i = 0; i < max; i++) {
         if (!get_next_line(fqr)) DIE("Improper fastq file");
         char *buf = getStartBuffer(fqr->buf);
         if (buf[0] != '@') DIE("Improper fastq file");
         if (!get_next_line(fqr)) DIE("Improper fastq file");
         int seqLen = getLengthBuffer(fqr->buf);
         bases += seqLen - 1;
         if (!get_next_line(fqr)) DIE("Improper fastq file");
         if (buf[0] != '+') DIE("Improper fastq file");
         if (!get_next_line(fqr)) DIE("Improper fastq file");
         if(seqLen != getLengthBuffer(fqr->buf)) DIE("Improper fastq file");
         if (fqr->fpos >= fqr->end_read) break;
    }
    int64_t fileSize = fqr->size;
    int64_t bytesRead = fqr->fpos - startPos;
    *estimatedReads = (fileSize * i + bytesRead - 1) / bytesRead;
    *estimatedBases = (fileSize * bases + bytesRead - 1) / bytesRead;
    destroy_fq_reader(fqr);
    return fileSize;
}
