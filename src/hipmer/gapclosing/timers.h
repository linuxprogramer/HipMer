#ifndef __TIMERS_H
#define __TIMERS_H

// for gaps.c
#define T_PROCESS_SRF 0
#define T_SRF_COUNTS 1
#define T_SRF_READS 2
#define T_CALC_SCAFF_DEPTH 3

#define T_PROCESS_CONTIGS 4

#define T_PROCESS_MERALIGNER 5
#define T_DIRECT_ALNS 6
#define T_PROJECTED_ALNS 7
#define T_GET_ALN 8
#define T_PUT_READ_ORIENT 9

#define T_PROCESS_FASTQ 10
#define T_FASTQ_FILE_READS 11
#define T_GET_READ_ORIENT 12
#define T_FASTQ_ATOMICS 13
#define T_GAP_COPY 14

#define T_CHECK_DHASH 15
#define T_PRINT_GAPS 16

#define T_BALANCE_GAPS 17

// for merauder4.c
#define T_CLOSE_GAPS 18
#define T_SPLINTING 19
#define T_MER_WALKS 20
#define T_PATCHING 21
#define T_SCAFF_TO_FASTA 22
#define T_OVERALL 23
#define T_CALC_ASSEMBLY_STATS 24

#define T_BARRIER 25

#define T_MMAP 26

#define UPC_TIMED_BARRIER do {                        \
        START_TIMER(T_BARRIER);                       \
        upc_barrier;                                  \
        stop_timer(T_BARRIER);                        \
    } while (0);

void init_timers(void);
#define START_TIMER(t) start_timer(t, #t)
void start_timer(int t, const char *desc);
void stop_timer(int t);
double get_elapsed_time(int t);
double get_timer_all_max(int t);
double get_timer_all_min(int t);
double get_timer_all_av(int t);
void print_timers(int serial_only);



#endif
