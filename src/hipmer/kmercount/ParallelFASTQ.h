#ifndef _PARALLEL_FASTQ_H_
#define _PARALLEL_FASTQ_H_

#include <mpi.h>
#include <string>
#include <vector>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <string.h>
#include <sstream>
#include "MPIUtils.h"

#include "../common/Buffer.h"
#include "../fqreader/fq_reader.h"

using namespace std;

class ParallelFASTQ
{
public:
    
    ParallelFASTQ(void)
        :nrecords(0), elapsed_t(0)
    {
        fqr = create_fq_reader();
        if (!fqr)
            DIE("Couldn't create fastq reader\n");
    }

    ~ParallelFASTQ()
    {
        close_fq(fqr);
        destroy_fq_reader(fqr);
    }

    void open(const char *filename, bool from_shm, long knownSize = -1)
    {
        if(fqr->f) close_fq(fqr);
        open_fq(fqr, filename, from_shm);
        assert(fqr->f);
    }
    
    size_t fill_block(vector<string> & seqs, vector<string> & quals, size_t count)
    {
        double t = MPI_Wtime();
        Buffer id = initBuffer(MAX_READ_NAME_LEN);
        Buffer seq = initBuffer(MAX_READ_LEN);
        Buffer qual = initBuffer(MAX_READ_LEN);;
        size_t records_read = 0;
        seqs.clear();
        quals.clear();
        for (int i = 0; i < count; i++) {
            if (!get_next_fq_record(fqr, id, seq, qual)) {
                elapsed_t += (MPI_Wtime() - t);
                freeBuffer(id);
                freeBuffer(seq);
                freeBuffer(qual);
                return records_read;
            }
            seqs.push_back(getStartBuffer(seq));
            quals.push_back(getStartBuffer(qual));
            records_read++;
            nrecords++;
        }
        elapsed_t += (MPI_Wtime() - t);
        freeBuffer(id);
        freeBuffer(seq);
        freeBuffer(qual);
        return records_read;
    }

    int64_t getTotalRecordsRead() { return nrecords; }
    double get_elapsed_time() { return elapsed_t; }

private:
    fq_reader_t fqr;
    int64_t nrecords;
    double elapsed_t;
};

#endif
