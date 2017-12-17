#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>

#include "Friends.h"
#include "KmerMiddle.hpp"
#include "readufx.h"
#ifdef DEBUGMMAP
	#include <mpi.h>
#endif
using namespace std;


#ifndef KMER_LENGTH
        #define KMER_LENGTH 51
#endif


FILE * OpenDebugFile(char * prefix, FILE * pFile, int my_rank)
{
    stringstream ss;
    string rank;
    ss << my_rank;
    ss >> rank;
    string ofilename = prefix;
    ofilename += rank;
    pFile  = fopen(ofilename.c_str(), "w");
    return pFile;
}

bool use_shm = false;
FILE * fptr;	// used when use_shm is false
char * addr; // used when use_shm is true
int64_t * sizeaddr;	// used when use_shm is true
int64_t remainingkmers; // used when use_shm is true
int64_t offset; 


// input parameters: filename, dsize, nprocs, my_rank
// output: myshare, dsize
FILE * UFXInitOpen(char * filename, int * dsize, int64_t *myshare, int nprocs, int my_rank, int64_t *nEntries, int from_shm)
{
	if (from_shm)
		use_shm = true;

    KmerMiddle::set_k(KMER_LENGTH);
    *dsize = sizeof(ufxpack);
    
    struct stat stufx;
    struct stat stsize;

	int localfsize = 0;

	if (use_shm) {
		stringstream ss;
		string rank;
		ss << my_rank;
		ss >> rank;
		string ofilename(filename);

		unsigned found = ofilename.find_last_of("/");
		ofilename = ofilename.substr(found+1);
		string fullname = ofilename + rank;
		string sizefile = fullname + string(".entries");	// filename<pid>.entries

		/* Read the size file first */
		int fs = shm_open(sizefile.c_str(), O_RDONLY , 0);
		fstat(fs, &stsize);
		if (fs <= 0) {
			printf("Error, shm_open failed %d %s when trying to open %s\n", errno, strerror(errno), sizefile.c_str());
			exit(-1);
		}
		sizeaddr = (int64_t*) mmap(NULL, stsize.st_size, PROT_READ, MAP_SHARED, fs, 0);
		if (sizeaddr == MAP_FAILED) {
			cout << "mmap open for read failed with " << errno << " - " << strerror(errno) << endl;
		}
		else {
			if(my_rank == 0) cout << "Opened mmap file " << sizefile << " that lists total ufx entries " << (*sizeaddr) << endl;
		}
		(*nEntries) =  (*sizeaddr);	// assign values
		close(fs);

		/* Now read the actual UFX file */	
		int fd = shm_open(fullname.c_str(), O_RDONLY , 0);
		fstat(fd, &stufx);
		if (fd <= 0) {
			printf("Error, shm_open failed %d %s when trying to open %s\n", errno, strerror(errno), fullname.c_str());
			exit(-1);
		}
		localfsize = stufx.st_size;
		addr = (char*) mmap(NULL , localfsize, PROT_READ , MAP_SHARED, fd, 0);
		if (addr == MAP_FAILED) {
			cout << "mmap open for read failed with " << errno << " - " << strerror(errno) << endl;
		}
		else {
			if(my_rank == 0) cout << "Opened mmap file " << fullname << " ufx binary" << endl;
		}
		close(fd);
	} else {
		stat(filename, &stufx);
	}
	
    int64_t filesize;
	if (use_shm) {
		*myshare  = localfsize / static_cast<int64_t>(*dsize );
		remainingkmers = *myshare;
		//int64_t offset = 0;

        #ifdef DEBUGMMAP
		// NB: this will fail if used with the current upc build, which does not link with mpi
		
		// When '-uses-mpi' is used to link a UPC application, the 'MPI_Init()' and 'MPI_Finalize()' calls are handled by the 
		// UPC runtime--these calls should not appear in client code. 
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Datatype datatype;
        MPI_Type_contiguous(sizeof(ufxpack), MPI_CHAR, &datatype );
        MPI_Type_commit(&datatype);

		// MPI_EXSCAN is used to perform a prefix reduction on data distributed across the group. 
		// The value in recvbuf on the process with rank 0 is undefined, and recvbuf is not signficant on process 0. 
		// The value in recvbuf on the process with rank 1 is defined as the value in sendbuf on the process with rank 0. 
		// For processes with rank i > 1, the operation returns, in the receive buffer of the process with rank i, 
		// the reduction of the values in the send buffers of processes with ranks 0,...,i-1 (inclusive)
        int64_t lengthuntil;
        MPI_Exscan(myshare, &lengthuntil, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
		if(my_rank == 0)	lengthuntil  = 0;

		MPI_File thefile;
		MPI_File_open(MPI_COMM_WORLD, "debug_ufx.bin", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);    

        MPI_File_set_view(thefile, lengthuntil * (*dsize), datatype, datatype, "external32", MPI_INFO_NULL);
        MPI_Status status;
		ufxpack * packed = (ufxpack *) addr;
        MPI_File_write_all(thefile, packed, (*myshare), datatype, &status);
        MPI_File_close(&thefile);
		MPI_Barrier(MPI_COMM_WORLD);
        #endif
		return NULL;
	} else {
		filesize = stufx.st_size;	
		if(my_rank  == 0) cout << "Filesize is " << filesize << " bytes" << endl;
		int64_t numentries = filesize / static_cast<int64_t>(*dsize );
		(*nEntries) = numentries;
		if(my_rank == 0) cout << "Number of records is " << numentries << endl;
		fptr = fopen(filename, "r");

		int64_t perproc = numentries / nprocs;
		int64_t begin = perproc * my_rank;
		if(my_rank == nprocs-1)
			*myshare = numentries - (nprocs-1)* (perproc);
		else		
			*myshare = perproc;

		fseek (fptr, begin * static_cast<int64_t>(*dsize ), SEEK_SET );
		return fptr;
	}
}


// inputs: f, dsize, requestedkmers, dmin
// outputs: kmersarr, counts, lefts, rights
// returns: number of k-mers read (can be less than requestedkmers if end of file)
int64_t UFXRead(int dsize, char *** kmersarr, int ** counts, char ** lefts, char ** rights, int64_t requestedkmers, int dmin, int reuse, int my_rank)
{
   	int64_t i;
	ufxpack *upack = NULL;
	int64_t totread = 0;

	if (!use_shm) {
		if(!fptr){
			cerr << "Thread " << my_rank << ": Problem reading binary input file\n";
			return 1;
		}
		upack = new ufxpack[requestedkmers];
		totread = fread(upack, dsize, requestedkmers, fptr);
	} else {
		upack = (ufxpack *) (addr + offset);	
        totread = std::min(remainingkmers, requestedkmers);
		remainingkmers -= totread; 
        offset = offset + totread * static_cast<int64_t>(dsize);
	}
    
    if(!reuse)  // OK in the last iteration too because the invariant (totread <= requestedkmers) holds
    {
        // (*kmersarr) is of type char**
        if (my_rank == 0) {
            cerr << "Thread 0: Allocating memory for " << totread << " reads: " << (totread * (sizeof(char*) + KMER_LENGTH+1 + sizeof(int) + 2) ) << " bytes" << endl;
        }
        (*kmersarr) = (char**) malloc(sizeof(char*) * totread);
        if (*kmersarr == NULL) { cerr << "Thread " << my_rank << ": Could not allocate memory for " << totread << " reads" << endl; exit(1); }
        for (i = 0; i < totread; i++) {
            (*kmersarr)[i] = (char*) malloc((KMER_LENGTH+1) * sizeof(char)); // extra character for NULL termination
            if ((*kmersarr)[i] == NULL) { cerr << "Thread " << my_rank << ": Could not allocate memory for read " << i << endl; exit(1); }
        }
    
        *counts = (int*) malloc(sizeof(int) * totread);
        *lefts = (char*) malloc(sizeof(char) * totread);
        *rights = (char*) malloc(sizeof(char) * totread);
        if (*counts == NULL || *lefts == NULL || *rights == NULL) { cerr << "Thread " << my_rank << ": Could not allocate memory for 3 * " << totread << " reads" << endl; exit(1); }
    }
    
	for(i=0; i< totread; ++i)
	{
		KmerMiddle kmer(upack[i].arr);
        
        // from C++11 standard:
        // 21.4.7.1 says that the pointer returned by c_str() must point to a buffer of length size()+1.
        // 21.4.5 says that the last element of this buffer must have a value of charT() -- in other words, the null character
		std::strcpy ((*kmersarr)[i], kmer.toString().c_str()); // (*kmersarr)[i] is of type char*
		(*counts)[i] = upack[i].count;
		if(upack[i].leftmax < dmin)
			(*lefts)[i] = 'X';
		else if(upack[i].leftmin < dmin)	
			(*lefts)[i] = upack[i].left;
		else	// both dmin < leftmin < leftmax
			(*lefts)[i] = 'F';

		if(upack[i].rightmax < dmin)
			(*rights)[i] = 'X';
		else if(upack[i].rightmin < dmin)	
			(*rights)[i] = upack[i].right;
		else	// both dmin < rightmin < rightmax
			(*rights)[i] = 'F';
        
        // fprintf(stderr, "K-mer is named (%c) %s (%c) with count %d\n", (*lefts)[i], (*kmersarr)[i], (*rights)[i] , (*counts)[i]);
	}
	if (!use_shm)
    	delete [] upack;

	return totread;
}

void DeAllocateAll(char *** kmersarr, int ** counts, char ** lefts, char ** rights, int64_t initialread)
{
   int64_t i;
    for (i = 0; i < initialread; i++)
        free((*kmersarr)[i]);
    free(*kmersarr);
    free(*counts);
    free(*lefts);
    free(*rights);
}

