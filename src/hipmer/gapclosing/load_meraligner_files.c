#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <upc.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <libgen.h>
#include <unistd.h>

#include "timers.h"
#include "utils.h"

int main(int argc, char** argv)
{
    if (argc != 3)
        return 0;    
    init_timers();
    char *dirname = strdup(argv[1]);
    int extension = atoi(argv[2]);
    char fname[PATH_MAX];
    sprintf(fname, "merAlignerOutput_%d_ReadFile%d", MYTHREAD, extension);
    char shm_fname[PATH_MAX];
    char disk_fname[PATH_MAX];
    sprintf(shm_fname, "/dev/shm/%s", fname);
    sprintf(disk_fname, "%s/%s", dirname, fname);
    serial_printf("Thread %d: Loading %s into %s\n", MYTHREAD, disk_fname, shm_fname);
    //fprintf(stderr, "Thread %d: Loading %s into %s\n", MYTHREAD, disk_fname, shm_fname);
    long file_size = get_file_size(disk_fname);
    FILE *shm_f = fopen_chk(shm_fname, "w+");
    FILE *f = fopen_chk(disk_fname, "r");
    START_TIMER(T_OVERALL);
    // read in 10 blocks so we can get timing
    long sub_block = INT_CEIL(file_size, 10);
    char *buf = malloc(sub_block);
    size_t bytes_read = 0;
    size_t bytes_written = 0;
    for (int i = 0; i < 10; i++) {
        long bytes = fread(buf, 1, sub_block, f);
        bytes_read += bytes; 
        bytes_written += fwrite(buf, 1, bytes, shm_f);
        serial_printf("%d ", i + 1);
    }
    fclose(f);
    fclose(shm_f);
    free(buf);
    if (bytes_read != bytes_written)
        WARN("Bytes read != bytes written: %ld ! %ld\n", bytes_read, bytes_written);
    serial_printf("\n");
    double gb_read = (double)bytes_read / ONE_GB;
    stop_timer(T_OVERALL);
    serial_printf("Read %.4f GB, %.4f GB/s\n", gb_read, gb_read / get_elapsed_time(T_OVERALL));
	upc_barrier;
    print_timers(1);
    return 0;
}
