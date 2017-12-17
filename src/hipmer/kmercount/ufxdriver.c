#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "readufx.h"

#if defined(USEMPI)
#include <mpi.h>
#elif defined(USEUPC)
#include <upc.h>
#endif

#define BLOCK2READ  1000000 // number of k-mers to read at each iteration by each processor

int main(int argc, char* argv[])
{
    int nprocs, myrank; 
#if defined(USEMPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#elif defined(USEUPC)
    nprocs = THREADS;
    myrank = MYTHREAD;
#else   // serial
    nprocs = 1;
    myrank = 0;
#endif
    
    if(argc < 3)
    {
        fprintf(stderr, "Usage: ./ufxdriver <ufxfilename> <dmin> <shm>\n");
    }
#ifdef DEBUG
    fprintf(stdout, "size_t is of size %ld\n", sizeof(size_t));
#endif
    
	int use_shm = 0;
	if (argc == 4 && strcmp(argv[3], "shm") == 0)
		use_shm = 1;
    // offsetting/fseeking to the file pointer done inside UFXInitOpen
    int64_t myshare; // how much data each process will "receive" in total
    int dsize;  // how many bytes each struct take
    int64_t total;	// total number of UFX entries in file

    UFXInitOpen(argv[1], &dsize, &myshare, nprocs, myrank, &total, use_shm);
#ifdef DEBUG
    fprintf(stderr, "UFX file opened, dsize: %d, myshare: %ld out of total %ld\n", dsize, myshare, total);
#endif
    
    char ** kmersarr;   // string array
    int * counts;
    char * lefts;
    char * rights;
	
    FILE * pFile;
    pFile = OpenDebugFile("kmerdebug", pFile, myrank);
    int realblock = (BLOCK2READ < myshare)? BLOCK2READ: myshare;

    int64_t readsofar = 0;
    int64_t initialread = UFXRead(dsize, &kmersarr, &counts, &lefts, &rights, realblock, atoi(argv[2]), 0, myrank);
    readsofar += initialread;
    for(int i=0; i< initialread; ++i)
        fprintf(pFile, "%c[%s]%c %d\n", lefts[i], kmersarr[i], rights[i], counts[i]);
    while( readsofar < myshare )
    {
        realblock = (BLOCK2READ < (myshare - readsofar))? BLOCK2READ: (myshare - readsofar);
        int64_t read = UFXRead(dsize, &kmersarr, &counts, &lefts, &rights, realblock, atoi(argv[2]), 1, myrank);
        readsofar += read;
        for(int i=0; i< read; ++i)
            fprintf(pFile, "%c[%s]%c %d\n", lefts[i], kmersarr[i], rights[i], counts[i]);
    }
    fclose(pFile);
    DeAllocateAll(&kmersarr, &counts, &lefts, &rights, initialread);
#if defined(USEMPI)
    MPI_Finalize();
#endif
    return 0;
}
