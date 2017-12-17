/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#include <string.h>
#include <upc.h>
#include <unistd.h>

#include "utils.h"

void get_mem_usage(memory_usage_t *mem)
{
    FILE *f = fopen_chk("/proc/self/statm", "r");
    if (fscanf(f, "%lu %lu %lu", &(mem->vm_size), &(mem->vm_rss), &(mem->share)) != 3) {
        WARN("could not read memory usage\n");
        memset(mem, 0, sizeof(memory_usage_t));
        return;
    }
    fclose(f);
    int ps = getpagesize();
    mem->vm_size *= ps;
    mem->vm_rss *= ps;
    mem->share *= ps;
}

int get_max_mem_usage_mb(void)
{
    int usage;
    FILE *f = fopen_chk("/proc/self/status", "r");
    char buf[1000];
    char units[5];
    while (fgets(buf, 999, f)) {
        if (strncmp(buf, "VmHWM:", 6) == 0) {
            sscanf(buf + 6, "%d %s\n", &usage, units);
            break;
        }
    }
    if (strcmp(units, "kB") == 0)
        usage /= 1024;
    return usage;
}

int bankers_round(double x)
{
    int rx = x;
    if (x - rx > 0.5)
        rx++;
    if (x - rx == 0.5 && rx % 2)
        rx++;
    return rx;
}

// reverse in place
char *reverse(char *s)
{
    int len = strlen(s);
    for (int i = 0; i < len / 2; i++) {
        char tmp = s[len - i - 1];
        s[len - i - 1] = s[i];
        s[i] = tmp;
    }
    return s;
}

int endswith(char *s1, char *s2)
{
    int s1l = strlen(s1);
    int s2l = strlen(s2);
    if (s1l < s2l) 
        return -1;
    return strcmp(s1 + strlen(s1) - strlen(s2), s2);
}

void switch_code(char *seq) {
    char *p = seq;
    while (*p) {
        switch (*p) {
        case 'a': *p = 't'; break;
        case 'c': *p = 'g'; break;
        case 'g': *p = 'c'; break;
        case 't': *p = 'a'; break;
        case 'A': *p = 'T'; break;
        case 'C': *p = 'G'; break;
        case 'G': *p = 'C'; break;
        case 'T': *p = 'A'; break;
        }
        p++;
    }
}

int get_tokens(char *s, char *tokens[], int *tot_len, int ntoks)
{
    char *s_begin = s;
    char *s_end;
    int tok_i = 0;
    *tot_len = 0;
    while ((s_end = strchr(s_begin, '\t')) != NULL) {
        s_end[0] = '\0';
        tokens[tok_i] = s_begin;
        tok_i++;
        *tot_len = s_end - s;
        if (tok_i == ntoks)
            break;
        s_begin = s_end + 1;
    }
    // get end of line?
    if (tok_i < ntoks && (s_end = strchr(s_begin, '\n')) != NULL) {
        s_end[0] = '\0';
        tokens[tok_i] = s_begin;
        *tot_len = s_end - s;
        tok_i++;
    }
    return tok_i;
}

static shared long _all_long[THREADS];

long reduce_long(long myval, upc_op_t op) {
    shared long *per_thread_long = upc_all_alloc(THREADS, sizeof(long));
    per_thread_long[MYTHREAD] = myval;
    bupc_all_reduce_allL(_all_long, per_thread_long, op, THREADS, 1, NULL, UPC_IN_MYSYNC|UPC_OUT_MYSYNC);
    upc_all_free(per_thread_long);
    return _all_long[MYTHREAD];
}

static shared int _all_int[THREADS];

int reduce_int(int myval, upc_op_t op) {
    shared int *per_thread_int = upc_all_alloc(THREADS, sizeof(int));
    per_thread_int[MYTHREAD] = myval;
    bupc_all_reduce_allI(_all_int, per_thread_int, op, THREADS, 1, NULL, UPC_IN_MYSYNC|UPC_OUT_MYSYNC);
    upc_all_free(per_thread_int);
	return _all_int[MYTHREAD];
}


/*
int pin_thread(pid_t pid, int cid) 
{
    cpu_set_t cpu_set;
    CPU_ZERO(&cpu_set);
    CPU_SET(cid, &cpu_set);
    if (sched_setaffinity(pid, sizeof(cpu_set), &cpu_set) == -1) {
        if (errno == 3) 
            WARN("%s, pid: %d", strerror(errno), pid);
        return -1;
    }
    return 0;
}
*/


