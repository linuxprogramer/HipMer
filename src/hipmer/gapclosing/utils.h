/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#ifndef __UTILS_H
#define __UTILS_H

#include <sys/time.h>
#include <stdio.h>
#include <errno.h>
#include <upc_collective.h>
#include <upc_tick.h>

#include "../common/common.h"

typedef struct {
    size_t vm_size;
    size_t vm_rss;
    size_t share;
} memory_usage_t;

#define CHECK_LINE_LEN(buf, buf_size, fname)                            \
    do {                                                                \
        if (strlen(buf) == buf_size - 1)                                \
            DIE("line too long in file %s: > %d", fname, buf_size - 1); \
    } while (0)

#define ONE_KB 1024
#define ONE_MB 1048576
#define ONE_GB 1073741824L

void get_mem_usage(memory_usage_t *mem);
int get_max_mem_usage_mb(void);

int bankers_round(double x);

char *reverse(char *s);
int endswith(char *s1, char *s2);
void switch_code(char *seq);

int get_tokens(char *s, char *tokens[], int *tot_len, int ntoks);

long reduce_long(long myval, upc_op_t op);
int reduce_int(int myval, upc_op_t op);

#endif
