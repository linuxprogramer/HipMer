/* UPC version of gap-closing, Steven Hofmeyr (shofmeyr@lbl.gov), Nov 2014.
 */

#ifndef __TRACING_H
#define __TRACING_H

#include <sys/time.h>
#include <stdio.h>
#include <errno.h>
#include <upc_collective.h>
#include <upc_tick.h>

#define tprintf_flush(fmt,...)                  \
    do {                                        \
        tprintf(fmt, ##__VA_ARGS__);            \
        tflush();                               \
    } while (0);

void tprintf(const char *fmt, ...);
void tflush(void);
void dbg(const char *fmt, ...);
void set_output_dir(const char *dirname);
void set_log_prefix(const char *logname);

#endif
