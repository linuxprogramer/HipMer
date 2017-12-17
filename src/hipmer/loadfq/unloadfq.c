#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glob.h>
#include <upc.h>

#include "../fqreader/fq_reader.h"

#define MAX_SEQ_FILES 12

int main(int argc, char** argv)
{
    if (argc == 1)
        return 0;    
    char *fname_glob = strdup(argv[1]);
    glob_t globbuf = {.gl_offs = MAX_SEQ_FILES};
    int err = glob(fname_glob, 0, 0, &globbuf);
    if (err)
        DIE("Could not obtain FASTQ file list %s: %s\n", fname_glob, strerror(err));
    if (!globbuf.gl_pathc)
        DIE("Couldn't find any matching files for %s\n", fname_glob);
    for (int fi = 0; fi < globbuf.gl_pathc; fi++) {
        char *fname = globbuf.gl_pathv[fi];
        serial_printf("Unloading %s from memory for %d threads\n", fname, THREADS);
        int err = unload_fq(fname);
        if (err < 0)
            serial_printf("Couldn't unload %s: %s\n", fname, strerror(-err));
    }
    return 0;
}
