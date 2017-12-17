#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Buffer.h"

// returns a new Buffer
Buffer initBuffer(size_t initSize) {
    if (initSize <= 0) initSize = 256;
    Buffer b = (_Buffer *) malloc(sizeof(_Buffer));
    if (b == NULL) { fprintf(stderr, "Could not allocate new Buffer!\n"); fflush(stderr); assert(0); exit(1); }
    b->buf = (char*) malloc(initSize);
    if (b->buf == NULL) { fprintf(stderr, "Could not allocate %ld bytes into Buffer!\n", initSize); fflush(stderr); assert(0); exit(1); }
    b->pos = b->len = 0;
    b->size = initSize;
    assert(b->len < b->size);
    b->buf[b->len] = '\0';
    assert(isValidBuffer(b));
    return b;
}
// returns the remaining size of the Buffer
size_t growBuffer(Buffer b, size_t appendSize) {
    assert(isValidBuffer(b));
    assert(b->size > 0);
    size_t requiredSize = b->len + appendSize + 1;
    if (requiredSize >= b->size) {
        while (requiredSize > b->size) {
            b->size *= 2;
        }
        assert(b->size >= requiredSize);
        b->buf = (char*) realloc(b->buf, b->size);
        if (b->buf == NULL)  { fprintf(stderr, "Could not reallocate %ld bytes into Buffer!", b->size); fflush(stderr); assert(0); exit(1); }
        b->buf[b->len] = '\0';
    }
    assert(b->len + appendSize < b->size);
    return b->size - b->len;
}

// grows the buffer to the maximum available size up to requestedSize
size_t growBufferMax(Buffer b, size_t requestedSize) {
    assert(isValidBuffer(b));
    assert(b->size > 0);
    size_t requiredSize = b->len + requestedSize + 1;
    char *newBuf = NULL;
    while ((newBuf = realloc(b->buf, requiredSize)) == NULL) {
        requiredSize = 3 * requiredSize / 4;
        if (requiredSize < b->size) { fprintf(stderr, "Could not growBufferMax to %ld (beyond %ld)\n", requestedSize, b->size); fflush(stderr); assert(0); exit(1); }
    }
    b->buf = newBuf;
    b->size = requiredSize;
    return b->size - b->len;
}

// sets the buffer len to 0, does not affect size
void resetBuffer1(Buffer b) {
    assert(b != NULL);
    assert(b->buf != NULL);
    assert(b->size > 0);
    b->pos = b->len = 0;
}
void resetBuffer(Buffer b) {
    resetBuffer1(b);
    assert(b->len == 0);
    assert(b->pos == 0);
    b->buf[b->len] = '\0';
}
char * resetRawBuffer(Buffer b, size_t newSize) {
    assert(isValidBuffer(b));
    resetBuffer(b);
    growBuffer(b, newSize);
    b->len = newSize;
    b->buf[b->len] = '\0';
    if (b->pos > b->len) b->pos = b->len;
    return b->buf;
}

// destroys a Buffer
void freeBuffer(Buffer b) {
    if (b == NULL) return;
    if (b->buf != NULL) free(b->buf);
    b->buf = NULL;
    b->pos = b->len = b->size = 0;
    free(b);
}

// releases the buffer memory (free now expected by user not freeBuffer)
char * releaseBuffer(Buffer b) {
    char *buf = b->buf;
    b->buf = NULL;
    freeBuffer(b);
    return buf;
}

int isValidBuffer(Buffer b) {
     assert(b != NULL) ;
     return (b != NULL) & (b->buf != NULL) & (b->size > 0) & (b->len >= b->pos) & (b->size >= b->len);
}


BufferList initBufferList(size_t initSize) {
	BufferList bl;
	bl.head = NULL;
	bl.tail = NULL;
	extendBufferList(bl, initSize);
	return bl;
}

void freeBufferList(BufferList bl) {
	struct _BufferList *ptr = bl.head;
	while (ptr != NULL) {
		struct _BufferList *copy = ptr;
		Buffer b = copy->buffer;
		freeBuffer(b);
		ptr = copy->next;
		copy->buffer = NULL;
		copy->next = NULL;
		free(copy);
	}
	bl.head = NULL;
	bl.tail = NULL;
}
	
Buffer extendBufferList(BufferList bl, size_t initSize) {
	struct _BufferList *newnode = (struct _BufferList *) calloc(1, sizeof(struct _BufferList));
	newnode->buffer = initBuffer(initSize);
	newnode->next = NULL;
	if (bl.tail) {
		bl.tail->next = newnode;
		bl.tail = newnode;
	} else {
		assert(!bl.head);
		bl.head = bl.tail = newnode;
	}
	return getBuffer(bl);
}

Buffer getBuffer(BufferList bl) {
	assert(bl.tail);
	assert(bl.tail->buffer);
	return bl.tail->buffer;
}

size_t getPosBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->pos;
}
size_t getLengthBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->len;
}
size_t getSizeBuffer(Buffer b) {
    return b->size;
}
char * getStartBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->buf;
}
char *getEndBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->buf + b->len;
}
char *getCurBuffer(Buffer b) {
    assert(isValidBuffer(b));
    return b->buf + b->pos;
}

size_t writeBuffer(Buffer b, const void *data, size_t size) {
    assert(isValidBuffer(b));
    memcpyBuffer(b, data, size);
    return size;
}

size_t readBuffer(Buffer b, void *data, size_t size) {
    assert(isValidBuffer(b));
    assert(b->pos <= b->len);
    size_t maxReadSize = b->len - b->pos;
    if (size >= maxReadSize) size = maxReadSize;
    if (size) {
        char *x = memcpy(data, b->buf + b->pos, size);
        assert(x == (char*) data);
        b->pos += size;
    }
    return size;
}

// writes to a buffer, reallocating, if necessary
// returns the length of the write
size_t printfBuffer(Buffer b, const char *fmt, ...)
{
    assert(isValidBuffer(b));
    size_t requiredSize = -1;
    assert(b != NULL);
    {
        va_list args;

        va_start(args, fmt);
        requiredSize = vsnprintf(NULL, 0, fmt, args);
        va_end(args);

        assert(requiredSize > 0);
        growBuffer(b, requiredSize+1);
    }
    assert(b->size - b->len > 0);

    va_list args;
    va_start(args, fmt);
    size_t len = vsnprintf(getEndBuffer(b), b->size - b->len, fmt, args);
    va_end(args);

    assert(len == requiredSize);
    b->len += len;
    assert(b->len < b->size);
    return len;
}

// reads one line or end of file, regardless of the line length
char *fgetsBuffer(Buffer b, int size, FILE *stream)
{
    assert(isValidBuffer(b));
    size_t oldLen = b->len;
    assert(oldLen < b->size);
    do {
        growBuffer(b, size);
        char *pos = fgets(getEndBuffer(b), size, stream);
        if (pos) {
            assert(pos == getEndBuffer(b));
            size_t readlen = strlen(pos);
            assert(readlen < size);
            assert(readlen > 0);
            b->len += readlen;
            assert(b->len < b->size);
            assert(pos + readlen == getEndBuffer(b));
            assert(*(getEndBuffer(b)) == '\0');
            if (b->buf[b->len-1] == '\n') break;
        } else {
            // nothing read
            assert(b->len < b->size);
            b->buf[b->len] = '\0';
            break;
        }
    } while(!feof(stream));
    // return the pointer to the beginning of this line
    return (b->len > oldLen) ? (getStartBuffer(b) + oldLen) : NULL;
}

void *memcpyBuffer(Buffer b, const void *src, size_t len)
{
    assert(isValidBuffer(b));
    growBuffer(b, len+1);
    void * ret = memcpy(getEndBuffer(b), src, len);
    b->len += len;
    assert(b->size > b->len);
    return ret;
}

char *strcpyBuffer(Buffer b, const char *src)
{
    assert(isValidBuffer(b));
    assert(src != NULL);
    size_t len = strlen(src);
    return strncpyBuffer(b, src, len);
}

char *strncpyBuffer(Buffer b, const char *src, size_t len)
{
    assert(isValidBuffer(b));
    growBuffer(b, len + 1);
    char *ret = strcpy(getEndBuffer(b), src);
    b->len += len;
    assert(b->size > b->len);
    return ret;
}

/* code to get fmemopen working on MAC / BSD */
#if defined (__APPLE__) && defined(__MACH__)
#include <sys/mman.h>

static int readfn(void *handler, char *buf, int size) {
  Buffer b = (Buffer) handler;
  size_t available = b->len - b->pos;
  
  if (size > available) {
    size = available;
  }
  memcpy(buf, b->buf + b->pos, sizeof(char) * size);
  b->pos += size;
  
  return size;
}

static int writefn(void *handler, const char *buf, int size) {
  Buffer b = (Buffer) handler;
  size_t available = b->size - b->len;

  if (size > available) {
    size = available;
  }
  memcpy(b->buf + b->len, buf, sizeof(char) * size);
  b->len += size;

  return size;
}

static fpos_t seekfn(void *handler, fpos_t offset, int whence) {
  size_t pos;
  Buffer b = (Buffer) handler;

  switch (whence) {
    case SEEK_SET: {
      if (offset >= 0) {
        pos = (size_t)offset;
      } else {
        pos = 0;
      }
      break;
    }
    case SEEK_CUR: {
      if (offset >= 0 || (size_t)(-offset) <= b->pos) {
        pos = b->pos + (size_t)offset;
      } else {
        pos = 0;
      }
      break;
    }
    case SEEK_END: pos = b->len + (size_t)offset; break;
    default: return -1;
  }

  if (pos > b->len) {
    return -1;
  }

  b->pos = pos;
  return (fpos_t)pos;
}

static int closefn(void *handler) {
  return 0; // freeing the Buffer should happen outside this scope, just as allocation 
}

FILE *fmemopen(void *buf, size_t size, const char *mode) {
  Buffer b = (Buffer) buf;

  // initialization and allocation of the buffer should have already happend
  assert(isValidBuffer(b));
  assert(b->size > 0);

  int reading = 0, writing = 0;
  int modlen = strlen(mode);
  assert(b->size == size);
    if (mode[0] == 'a') {
        writing = 1;
        if (strlen(mode) > 1 && mode[1] == '+') {
            reading = 1;
        }
    } else {
        if (mode[0] == 'r') {
            /* reading */
            reading = 1;
        } else if (mode[0] == 'w') {
            writing = 1;
        if (mode[modelen-1] == '+') {
            reading = 1;
            writing = 1;
        }
    }


  // funopen's man page: https://developer.apple.com/library/mac/#documentation/Darwin/Reference/ManPages/man3/funopen.3.html
  return funopen(b, reading?readfn:NULL, writing?writefn:NULL, seekfn, closefn);
}
#endif /* Helper fmemopen methods for Mac OS X */


/* returns a FILE object linked to the Buffer */
FILE *getBufferFile(Buffer b, const char *mode) {
    // initialization and allocation of the buffer should have already happend
    assert(isValidBuffer(b));
    assert(b->size > 0);
    
    /* reset pos (read position) or len (write position) for the buffer if not appending */
    int modelen = strlen(mode);
    if (mode[0] == 'r') {
        /* reading */
        b->pos = 0;
    } else if (mode[0] == 'w') {
        b->len = 0;
    }
    if (mode[modelen-1] == '+') {
        if (mode[0] == 'w' || mode[0] == 'a') {
            b->pos = 0;
        } else {
            assert(mode[0] == 'r');
            b->len = 0;
        }
    }
    
    FILE *f = fmemopen(b->buf, b->len, mode);
    setbuf(f, NULL);
    return f;
}

void setBufferForFile(Buffer b, FILE *f) {
    setbuffer(f, b->buf, b->size);
}

// removes a trailing new line, if it exists
// returns the possibly new length
size_t chompBuffer(Buffer b)
{
    assert(isValidBuffer(b));
    assert(b->buf[b->len] == '\0');
    if (b->buf[b->len - 1] == '\n') {
        b->buf[b->len - 1] = '\0';
        b->len--;
    }
    return b->len;
}

void swapBuffer(Buffer *a, Buffer *b) {
    Buffer tmp = *a;
    *a = *b;
    *b = tmp;
}
