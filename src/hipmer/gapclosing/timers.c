#include <string.h>

#include "utils.h"
#include "tracing.h"
#include "timers.h"

#define TICKS_TO_S(t) ((double)upc_ticks_to_ns(t) / 1000000000.0)
#define ELAPSED_TIME(t) TICKS_TO_S(upc_ticks_now() - (t))

typedef struct {
    char *description;
    upc_tick_t t_start;
    upc_tick_t t_stop;
    upc_tick_t t_elapsed;
} upc_timer_t;

#define MAX_TIMERS 200

static upc_timer_t _timers[MAX_TIMERS];

void init_timers(void)
{
    for (int i = 0; i < MAX_TIMERS; i++)
        memset(&_timers[i], 0, sizeof(upc_timer_t));
}

void start_timer(int t, const char *desc)
{
    _timers[t].t_start = upc_ticks_now();
    if (!_timers[t].description)
        _timers[t].description = strdup(desc);
#ifdef DEBUG
    tprintf("Starting timer for %s\n", desc);
#endif
}

void stop_timer(int t)
{
    _timers[t].t_stop = upc_ticks_now();
    double duration = (_timers[t].t_stop - _timers[t].t_start);
    _timers[t].t_elapsed += duration;
#ifdef DEBUG
    tprintf("Stopped timer for %s %0.3f (%0.3f)\n", _timers[t].description, TICKS_TO_S(duration), TICKS_TO_S(_timers[t].t_elapsed));
#endif
}

double get_elapsed_time(int t)
{ 
    return TICKS_TO_S(_timers[t].t_elapsed);
}

double get_timer_all_max(int t)
{
    return TICKS_TO_S(reduce_long(_timers[t].t_elapsed, UPC_MAX));
}

double get_timer_all_min(int t)
{
    return TICKS_TO_S(reduce_long(_timers[t].t_elapsed, UPC_MIN));
}

double get_timer_all_av(int t)
{
    return TICKS_TO_S(reduce_long(_timers[t].t_elapsed, UPC_ADD)) / THREADS;
}

void print_timers(int serial_only)
{
    //if (!serial_only)
    //    tprintf("%-27s %10s %10s %10s %10s %10s\n", "Timers", "thread", "min", "av", "max", "bln");
    serial_printf("  %-25s %10s %10s %10s %10s\n", "Timers", "min", "av", "max", "bln");
    for (int i = 0; i < MAX_TIMERS; i++) {
//        tprintf_flush("Timer: print_timers loop: %d\n");
        double t_min = get_timer_all_min(i);
        double t_av = get_timer_all_av(i);
        double t_max = get_timer_all_max(i);
        double bln = 0.0;
        if (t_max > 0.0) bln = t_av / t_max;
        if (!_timers[i].description)
            continue;
        serial_printf("  %-25s %10.2f %10.2f %10.2f %10.2f\n", _timers[i].description, 
                      t_min, t_av, t_max, bln);
        //if (!serial_only)
        //   tprintf_flush("  %-25s %10.2f %10.2f %10.2f %10.2f %10.2f\n", _timers[i].description, 
        //            TICKS_TO_S(_timers[i].t_elapsed), t_min, t_av, t_max, bln);
    }
}
