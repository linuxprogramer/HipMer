/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#include <sys/stat.h>
#include <string.h>
#include <upc.h>
#include <unistd.h>

#include "../common/common.h"
#include "tracing.h"
#include "utils.h"

static char _output_dir[MAX_FILE_PATH] = ".";

void set_output_dir(const char *dirname)
{
    strcpy(_output_dir, dirname);
}

static FILE *my_file = NULL;
static char LOG_PREFIX[MAX_FILE_PATH] = "";
#define DEFAULT_LOG_PREFIX "gapclosing"

void set_log_prefix(const char *name)
{
    strcpy(LOG_PREFIX, name);
    get_rank_path(LOG_PREFIX, MYTHREAD);
}


void tprintf(const char *fmt, ...)
{
    if (!my_file) {
        char buf[MAX_FILE_PATH];
        if (LOG_PREFIX[0] == '\0') {
            set_log_prefix(DEFAULT_LOG_PREFIX);
        }
        sprintf(buf, "%s/%s_%d.log", _output_dir, LOG_PREFIX, MYTHREAD);
        my_file = fopen_chk(buf, "w");
    }
    va_list args;
    va_start(args, fmt);
    vfprintf(my_file, fmt, args);
    va_end(args);
#ifdef DEBUG
    // always flush messages if debugging
    tflush();
#endif
}

void tflush(void)
{
    fflush(my_file);
}

static FILE *dbg_file = NULL;

void dbg(const char *fmt, ...)
{
    if (!dbg_file) {
        char buf[MAX_FILE_PATH];
        if (LOG_PREFIX[0] == '\0') {
            set_log_prefix(DEFAULT_LOG_PREFIX);
        }
        sprintf(buf, "%s/%s_%d.dbg", _output_dir, LOG_PREFIX, MYTHREAD);
        //sprintf(buf, "/dev/shm/gapclosing_%d.dbg", MYTHREAD);
        dbg_file = fopen_chk(buf, "w");
    }
    va_list args;
    va_start(args, fmt);
    vfprintf(dbg_file, fmt, args);
    va_end(args);
#ifdef DEBUG
    // always flush if debugging
    fflush(dbg_file);
#endif
}



