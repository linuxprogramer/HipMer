#include "HipMERdefines.h"
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <vector>
#include <sstream> 
#include <limits>
#include <array>
#include <map>
#include <tuple>
#if GCC_VERSION < 40300
	#include <tr1/random>
#else
	#include <random>
#endif
#include <mpi.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "MPIUtils.h"
#include "DataStructures/hyperloglog.hpp"
extern "C" {
#ifdef HIPMER_BLOOM64
#include "DataStructures/libbloom/bloom64.h"
#else
#include "DataStructures/libbloom/bloom.h"
#endif
}
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "Deleter.h"
#include "ParallelFASTQ.h"
#include "Friends.h"
#include "MPIType.h"
#include "SimpleCount.h"
#include "FriendsMPI.h"
#include "omp.h"
#include "Pack.h"
#include "../common/Buffer.h"
#ifdef KHASH
	#include "khash.hh"
#endif

#include <getopt.h>

#ifdef HEAVYHITTERS
#ifndef MAXHITTERS
#define MAXHITTERS 32000
#endif
SimpleCount<Kmer> heavyhitters(MAXHITTERS);
UFX2ReduceObj * Frequents;
size_t nfreqs;
#endif

using namespace std;

//#define BENCHMARKONLY // does not output a ufx file when defined
#define MEGA 1000000.0
#define MILLION 1000000
#define PROCESSRAM 128 * 1024 * 1024	// provisioning for 128 MB RAM per core for data I/O buffers
#define COUNT_THRESHOLD 3000000
#define COUNT_THRESHOLD_HIGH 33000000
#define HIGH_BIN 10
#define HIGH_NUM_BINS ((COUNT_THRESHOLD_HIGH-COUNT_THRESHOLD)/HIGH_BIN)

#ifndef MAX_ALLTOALL_MEM
#define MAX_ALLTOALL_MEM 134217728
#endif

int ERR_THRESHOLD;
typedef array<int,4> ACGT;
int nprocs;
int myrank;
int64_t nonerrorkmers;
int64_t kmersprocessed;
int64_t readsprocessed;
int totaltime;
typedef array<uint64_t, KMERLONGS> MERARR;


typedef tuple<ACGT,ACGT,int> KmerCountType;
#ifdef KHASH
typedef khmap_t< MERARR, KmerCountType > KmerCountsType;
#else
typedef map< MERARR, KmerCountType > KmerCountsType;
#endif
KmerCountsType kmercounts;

// Kmer is of length k
// HyperLogLog counting, bloom filtering, and std::maps use Kmer as their key
size_t ParseNPack(vector<string> & seqs, vector<string> & quals, vector< vector<Kmer> > & outgoing,
                vector<vector<array<char,2>>> & extquals, vector<vector<array<char,2>>> & extseqs, int pass, size_t offset)
{
	assert(pass == 2);
	MPI_Pcontrol(1,"ParseNPack");
	char special = phred_encoding+2;
	size_t nreads = seqs.size();
	size_t maxsending = 0;
	size_t bytesperkmer = Kmer::numBytes();
	size_t bytesperentry = bytesperkmer + 4;

	for(size_t i=offset; i< nreads; ++i)
	{
		size_t found;
#ifdef TRIMSPECIAL
		found = quals[i].find_first_of(special);
		if(found == string::npos)
			found = seqs[i].length();	// remember that the last valid position is length()-1
#else
	        found = seqs[i].length();
#endif
		// skip this sequence if the length is too short
		if (seqs[i].length() <= KMER_LENGTH) {
			//cerr << "seq is too short (" << seqs[i].length() << " < " << KMER_LENGTH << " : " << seqs[i] << endl;
			continue;
		}
		kmersprocessed += (seqs[i].length()-KMER_LENGTH+1);
        
		for(size_t j=0; j< seqs[i].length()-KMER_LENGTH+1; ++j)
		{
			size_t sending = PackEnds(seqs[i], quals[i], j, outgoing, extquals, extseqs, pass, found);
			if (sending > maxsending) maxsending = sending;
		}
		if (maxsending * bytesperentry >= MAX_ALLTOALL_MEM / nprocs) return i+1;
	}
	MPI_Pcontrol(-1,"ParseNPack");
	return nreads;
}

void Exchange(vector< vector<Kmer> > & outgoing, vector<vector<array<char,2>>> & extquals, vector<vector<array<char,2>>> & extseqs,
              vector<Kmer> & mykmers, vector<array<char,2>> & myquals, vector<array<char,2>> & myseqs, int pass, Buffer scratch1, Buffer scratch2)
{
	assert(pass == 2);
	MPI_Pcontrol(1,"Exchange");

	// serialize k-mer
	size_t bytesperkmer = Kmer::numBytes();
	size_t bytesperentry = bytesperkmer + 4;
	int * sendcnt = new int[nprocs];
	for(int i=0; i<nprocs; ++i) {
		sendcnt[i] = (int) outgoing[i].size() * bytesperentry;
		assert( outgoing[i].size() == extquals[i].size() );
		assert( outgoing[i].size() == extseqs[i].size() );
        }
	int * sdispls = new int[nprocs];
	int * rdispls = new int[nprocs];
	int * recvcnt = new int[nprocs];
	CHECK_MPI( MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD) );  // share the request counts
	sdispls[0] = 0;
	rdispls[0] = 0;
	for(int i=0; i<nprocs-1; ++i) {
		if (sendcnt[i] < 0 || recvcnt[i] < 0) {
			cerr << myrank << " detected overflow in Alltoall" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		sdispls[i+1] = sdispls[i] + sendcnt[i];
		rdispls[i+1] = rdispls[i] + recvcnt[i];
		if (sdispls[i+1] < 0 || rdispls[i+1] < 0) {
			cerr << myrank << " detected overflow in Alltoall" << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
	uint64_t totsend = accumulate(sendcnt,sendcnt+nprocs, static_cast<uint64_t>(0));
	uint64_t totrecv = accumulate(recvcnt,recvcnt+nprocs, static_cast<uint64_t>(0));

#ifdef DEBUG
	uint64_t sendmin, sendmax;
	CHECK_MPI( MPI_Reduce(&totsend, &sendmin, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&totsend, &sendmax, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD) );
	if(myrank == 0) cout << "Max send: " << sendmax << ", Min send: " << sendmin << " totsend: " << totsend << endl;
#endif

	growBuffer(scratch1, sizeof(uint8_t) * totsend);
	uint8_t * sendbuf = (uint8_t*) getStartBuffer(scratch1);
	for(int i=0; i<nprocs; ++i)  {
		size_t nkmers2send = outgoing[i].size();
		uint8_t * addrs2fill = sendbuf+sdispls[i];
		for(size_t j=0; j< nkmers2send; ++j) {
			assert(addrs2fill == sendbuf+sdispls[i] + j*bytesperentry);
			(outgoing[i][j]).copyDataInto( addrs2fill );
			char *ptr = ((char*) addrs2fill) + bytesperkmer;
			ptr[0] = extquals[i][j][0];
	                ptr[1] = extquals[i][j][1];
			ptr[2] = extseqs[i][j][0];
			ptr[3] = extseqs[i][j][1];
			addrs2fill += bytesperentry;
		}
		outgoing[i].clear();
		extquals[i].clear();
		extseqs[i].clear();
	}

	growBuffer(scratch2, sizeof(uint8_t) * totrecv);
	uint8_t * recvbuf = (uint8_t*) getStartBuffer(scratch2);

#ifdef DEBUG
	uint64_t globalmin, globalmax;
	CHECK_MPI( MPI_Reduce(&totrecv, &globalmin, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Reduce(&totrecv, &globalmax, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD) );
	if(myrank == 0) cout << "Max recv: " << globalmax << ", Min recv: " << globalmin << " totrecv: " << totrecv << endl;
	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
#endif
	CHECK_MPI( MPI_Alltoallv(sendbuf, sendcnt, sdispls, MPI_BYTE, recvbuf, recvcnt, rdispls, MPI_BYTE, MPI_COMM_WORLD) );

	uint64_t nkmersrecvd = totrecv / bytesperentry;
	for(uint64_t i= 0; i < nkmersrecvd; ++i) {
		Kmer kk;
		kk.copyDataFrom(recvbuf + (i * bytesperentry));	
		mykmers.push_back(kk);
		char *ptr = ((char*) recvbuf) + (i * bytesperentry) + bytesperkmer;
		array<char,2> qualexts = { ptr[0], ptr[1] };
		array<char,2> seqexts = { ptr[2], ptr[3] };
            	myquals.push_back(qualexts);
            	myseqs.push_back(seqexts);
	}

        DeleteAll(rdispls, sdispls, recvcnt, sendcnt);

	if(myrank == 0) cout << "exchanged " << totsend << ", " << totrecv << endl;
	MPI_Pcontrol(-1,"Exchange");
}


void InsertIntoHLL(vector<string> & seqs, vector<string> & quals, HyperLogLog & hll)
{
    char special = phred_encoding+2;
    size_t locreads = seqs.size();
    readsprocessed += locreads;
    
    MPI_Pcontrol(1,"HLL_Parse");
    for(size_t i=0; i< locreads; ++i)
    {
#ifdef TRIMSPECIAL
	size_t found = quals[i].find_first_of(special);
	if(found == string::npos)
        found = seqs[i].length();	// remember that the last valid position is length()-1
#else
        size_t found = seqs[i].length();
#endif
		
        if(found >= KMER_LENGTH) // otherwise size_t being unsigned will underflow
        {
            for(size_t j=0; j<= found-KMER_LENGTH; ++j)
            {
                string kmerstr = seqs[i].substr(j, KMER_LENGTH);
                for (auto & c: kmerstr) c = toupper(c);	// convert all to uppercase
                size_t Nfound=kmerstr.find('N');
                if (Nfound!=std::string::npos) continue;	// if there is an 'N', toss it
                Kmer mykmer(kmerstr.c_str());
                Kmer lexsmall =  mykmer.rep();
                kmerstr = lexsmall.toString();	// update string with its lexicographically smaller version
                hll.add(kmerstr.c_str(), kmerstr.size());
      	#ifdef HEAVYHITTERS 
                heavyhitters.Push(lexsmall);
	#endif
            }
        }
    }
    MPI_Pcontrol(-1,"HLL_Parse");
    seqs.clear();
    quals.clear();
}

void ProudlyParallelCardinalityEstimate(const vector<filedata> & allfiles, double & cardinality, bool use_shm)
{
	HyperLogLog hll(12);

	for(auto itr= allfiles.begin(); itr != allfiles.end(); itr++) {
	        ParallelFASTQ pfq;
		pfq.open(itr->filename, use_shm, itr->filesize);
		// The fastq file is read line by line, so the number of records processed in a block 
		// shouldn't make any difference, so we can just set this to some arbitrary value.
		// The value 262144 is for records with read lengths of about 100.
		size_t upperlimit = 262144;
		size_t totalsofar = 0;
		vector<string> seqs;
		vector<string> quals;
		vector<Kmer> mykmers;
		int iterations = 0;

		while (1) {
			MPI_Pcontrol(1,"FastqIO");
			size_t fill_status = pfq.fill_block(seqs, quals, upperlimit);
			// Sanity checks
			assert(seqs.size() == quals.size());
			for (int i = 0; i < seqs.size(); i++) {
				if (seqs[i].length() != quals[i].length()) {
					fprintf(stderr, "sequence length %ld != quals length %ld\n%s\n", 
							seqs[i].length(), quals[i].length(), seqs[i].c_str());
					MPI_Finalize();
					exit(1);
				}
			}
			MPI_Pcontrol(-1,"FastqIO");
			if (!fill_status)
				break;
			// this clears the vectors
			InsertIntoHLL(seqs, quals, hll);
		}
       
		double t = pfq.get_elapsed_time();
		double tot_t = 0;
		CHECK_MPI( MPI_Reduce(&t, &tot_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
		if (myrank == 0) {
			int num_ranks;
			double t2 = pfq.get_elapsed_time();
			CHECK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &num_ranks) );
			cout << "Total time taken for FASTQ reads is " << (tot_t / num_ranks) << ", elapsed " << t2 << endl;
        	}
	}

	// using MPI_UNSIGNED_CHAR because MPI_MAX is not allowed on MPI_BYTE
	int count = hll.M.size();
	CHECK_MPI( MPI_Allreduce(MPI_IN_PLACE, hll.M.data(), count, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Allreduce(MPI_IN_PLACE, &readsprocessed, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD) );

	cardinality = hll.estimate();
	if(myrank == 0) {
		cout << "Embarrassingly parallel k-mer count estimate is " << cardinality << endl;
		cout << "Total reads processed over all processors is " << readsprocessed << endl;
	}
#ifdef HEAVYHITTERS
	MPI_Pcontrol(1,"HeavyHitters");
	ParallelAllReduce(heavyhitters);
	nfreqs = heavyhitters.Size();
	heavyhitters.CreateIndex(); // make it indexible
	Frequents = new UFX2ReduceObj[nfreqs]();    // default initialize
	MPI_Pcontrol(-1,"HeavyHitters");
#endif

	MPI_Barrier(MPI_COMM_WORLD);
	cardinality /= static_cast<double>(nprocs);	// assume a balanced distribution
	cardinality *= 1.1;	// 10% benefit of doubt
}

class KmerInfo {
public:
    typedef array<char,2> TwoChar;
private:
    Kmer kmer;
    TwoChar quals,seqs;
public:
    KmerInfo() {}
    KmerInfo(Kmer k, TwoChar q, TwoChar s): kmer(k), quals(q), seqs(s) {}
    KmerInfo(const KmerInfo &copy) {
        kmer = copy.kmer;
        quals = copy.quals;
        seqs = copy.seqs;
    }
    int write(FILE *f) {
        return fwrite(this, sizeof(*this), 1, f);
        
    }
    int read(FILE *f) {
        return fread(this, sizeof(*this), 1, f);
    }
    // returns true if already exists, false if otherwise
    bool check(struct bloom *bm) {
        bool found = false;
        MPI_Pcontrol(1,"BloomFilter");

        string kmerstr = kmer.toString();
        if(bloom_check(bm, kmerstr.c_str(), kmerstr.size()) == 1)
        {
            MERARR karr = kmer.getArray();
            auto got = kmercounts.find (karr);  // kmercounts is a global variable
            if(got == kmercounts.end())
            {
                ACGT nulltuple = {0,0,0,0};
            #ifdef KHASH
                kmercounts.insert(karr, make_tuple(nulltuple,nulltuple,0));
            #else
                kmercounts.insert(make_pair(karr, make_tuple(nulltuple,nulltuple,0)));
            #endif
            }
            found = true;
        }
        else
        {
            bloom_add(bm, kmerstr.c_str(), kmerstr.size());
        }
        MPI_Pcontrol(-1,"BloomFilter");
        if (found)
            includeCount();
        return found;
    }

    bool includeCount() {
        MERARR karr = kmer.getArray();
        auto got = kmercounts.find (karr);  // kmercounts is a global variable
        return includeCount(got);
    }

    bool includeCount(KmerCountsType::iterator got) {
        MPI_Pcontrol(1,"HashTable");
        bool inserted = false;
        if(got != kmercounts.end()) // don't count anything else
        {
            if (quals[0] >= 0 && quals[1] >= 0) {
                // count the kmer in mercount
            #ifdef KHASH
                ++(get<2>(*got));  // ::value returns (const valtype_t &) but ::* returns (valtype_t &), which can be changed
            #else
                ++(get<2>(got->second)); // increment the counter regardless of quality extensions
            #endif
            } else {
                // ignore the mercount, but still use the exteions, if good enough
                quals[0] = 0 - quals[0];
                quals[1] = 0 - quals[1];
            }
            if(quals[0] >= phred_encoding+extensionCutoff)
            {
            #ifdef KHASH
                Increment1stBasedOnExt(*got, seqs[0]);
            #else
                Increment1stBasedOnExt(got->second, seqs[0]);
            #endif
            }
            if(quals[1] >= phred_encoding+extensionCutoff)
            {
            #ifdef KHASH
                Increment2ndBasedOnExt(*got, seqs[1]);
            #else
                Increment2ndBasedOnExt(got->second, seqs[1]);
            #endif
            }
            inserted = true;
        }
        MPI_Pcontrol(-1,"HashTable");
        return inserted;
    }
};

// at this point, no kmers include anything other than uppercase 'A/C/G/T'
void DealWithInMemoryData(vector<Kmer> & mykmers, int pass, struct bloom * bm, vector<array<char,2>> & myquals, vector<array<char,2>> & myseqs, FILE *possibleErrors)
{
    // store kmer & extensions with confirmaton of bloom
    // store to disk possible kmer with (first) insert into bloom
    assert(possibleErrors != NULL);
    assert(pass == 1);
    if(pass == 1)
    {
        MPI_Pcontrol(1,"BloomFilter");
        size_t count = mykmers.size();
        for(size_t i=0; i < count; ++i)
        {
            KmerInfo ki(mykmers[i], myquals[i], myseqs[i]);
            if (ki.check(bm)) {
                // found and included this kmer
            } else {
                // store possible kmer in case it turns out to be good
                ki.write(possibleErrors);
            }
        }
        MPI_Pcontrol(-1,"BloomFilter");
    }
}


void ProcessFiles(const vector<filedata> & allfiles, int pass, double & cardinality, bool use_shm)
{
    struct bloom * bm = NULL;
    FILE *possibleErrors = NULL;
    assert(pass == 1);
    
    char possibleErrorFile[MAX_FILE_PATH];
    char *pefilename = possibleErrorFile;
    if (use_shm) {
        sprintf(possibleErrorFile, "/dev/shm/possibleErrorKmers_%d", myrank);
    } else {
        sprintf(possibleErrorFile, "possibleErrorKmers_%d", myrank);
        pefilename = get_rank_path(possibleErrorFile, myrank);
    }
    possibleErrors = fopen_chk(pefilename, "w+"); // create, truncate then open for writing and reading
    
    Buffer scratch1 = initBuffer(MAX_ALLTOALL_MEM);
    Buffer scratch2 = initBuffer(MAX_ALLTOALL_MEM);
    Buffer peBuffer = initBuffer(0);
    if (!use_shm) {
        // allocate a large buffer for fileops
        growBufferMax(peBuffer, MAX_ALLTOALL_MEM);
        setBufferForFile(peBuffer, possibleErrors);
    }
    if(pass == 1)
    {
        unsigned int random_seed = 0xA57EC3B2;
        const double desired_probability_of_false_positive = 0.05;
        bm = (struct bloom*) malloc(sizeof(struct bloom));
		//初始化bloom filter的大小
        bloom_init(bm, cardinality, desired_probability_of_false_positive);
        
        if(myrank == 0)
        {
            cout << "Table size is: " << bm->bits << " bits" << endl;
            cout << "Optimal number of hash functions is : " << bm->hashes << endl;
        }

	// Only perform exchange in pass 1!!

	double pfqTime = 0.0;
	auto files_itr = allfiles.begin();
	while(files_itr != allfiles.end())
	{
            ParallelFASTQ *pfq = new ParallelFASTQ();
            pfq->open(files_itr->filename, use_shm, files_itr->filesize);
            files_itr++;
            // once again, arbitrarily chosen - see ProudlyParallelCardinalityEstimate
            size_t upperlimit = 262144;
	
            vector<string> seqs;
            vector<string> quals;
            vector< vector<Kmer> > outgoing(nprocs);
            vector< vector<array<char,2> > > extquals(nprocs);
            vector< vector<array<char,2> > > extseqs(nprocs);
            
            vector<Kmer> mykmers;
            vector<array<char,2>> myquals;
            vector<array<char,2>> myseqs;

            double t01 = MPI_Wtime();

            int moreflags[3], allmore2go[3], anymore2go;
            int &moreSeqs = moreflags[0], &moreToExchange = moreflags[1], &moreFiles = moreflags[2];
            int &allmoreSeqs = allmore2go[0], &allmoreToExchange = allmore2go[1], &allmoreFiles = allmore2go[2];
            moreSeqs = moreToExchange = 1;
            size_t fill_status;
            do
            {
			MPI_Pcontrol(1,"FastqIO");
			if (pfq && !moreSeqs)
			{
				assert(pfq != NULL);
            			double t = pfq->get_elapsed_time();
				pfqTime += t;
				delete pfq;
				pfq = NULL;
				if (files_itr != allfiles.end()) {
					pfq = new ParallelFASTQ();
					pfq->open(files_itr->filename, use_shm, files_itr->filesize);
					files_itr++;
					moreSeqs = 1;
				}
			}
			moreFiles = (files_itr != allfiles.end());
			if(moreSeqs)
			{
				fill_status = pfq->fill_block(seqs, quals, upperlimit);  // file_status is 0 if fpos >= end_fpos
			}
			moreSeqs = (fill_status > 0);
			MPI_Pcontrol(-1,"FastqIO");

			size_t offset = 0;
			do {
				offset = ParseNPack(seqs, quals, outgoing, extquals, extseqs, 2, offset);    // no-op if seqs.size() == 0
				if (offset == seqs.size()) {
					seqs.clear();	// no need to do the swap trick as we will reuse these buffers in the next iteration
					quals.clear();
					offset = 0;
				}
				Exchange(outgoing, extquals, extseqs, mykmers, myquals, myseqs, 2, scratch1, scratch2); // outgoing arrays will be all empty, shouldn't crush
				DealWithInMemoryData(mykmers, pass, bm, myquals, myseqs, possibleErrors);   // we might still receive data even if we didn't send any
				moreToExchange = offset < seqs.size();
				mykmers.clear();
				myquals.clear();
				myseqs.clear();

//				more2go = (fill_status > 0) | (offset != seqs.size()) | (files_itr != allfiles.end());    // this proc has still data to process
				CHECK_MPI( MPI_Allreduce(moreflags, allmore2go, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD) );
				if (myrank == 0) {
					cout << "active ranks moreSeqs: " << allmoreSeqs
						<< " moreToExchange: " << allmoreToExchange
						<< " moreFiles: " << allmoreFiles  
						<< ", rank " << myrank 
						<< " moreSeqs: " << moreSeqs 
						<< " moreToExchange: " << moreToExchange
						<< " moreFiles: " << moreFiles <<  endl;
				}
			} while (moreToExchange);
			anymore2go = allmoreSeqs + allmoreToExchange + allmoreFiles;
            } while(anymore2go);

            if (pfq) {
                double t = pfq->get_elapsed_time();
                pfqTime += t;
                delete pfq;
                pfq = NULL;
            }

            double tot_t = 0;
            CHECK_MPI( MPI_Reduce(&pfqTime, &tot_t, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) );
            if (myrank == 0) {
                int num_ranks;
                CHECK_MPI( MPI_Comm_size(MPI_COMM_WORLD, &num_ranks) );
                cout << "Total time taken for FASTQ reads is " << (tot_t / num_ranks) << ", elapsed " << pfqTime << endl;
            }
            double t02 = MPI_Wtime();
	
            if(myrank == 0) 
            {
			cout << "PASS "<< pass << endl;
			cout << "Read/distributed/processed reads of " << (files_itr == allfiles.end() ? " ALL files " : files_itr->filename) << " in " << t02-t01 << " seconds" << endl;
            }
        }	// end_of_loop_over_all_files
	
        bloom_free(bm);
        free(bm);
        bm = NULL;

	freeBuffer(scratch1);
	freeBuffer(scratch2);

        MPI_Pcontrol(1,"ReprocessFirstObservations");
        KmerInfo ki;

        // rewind possible Errors
        size_t peFileSize = ftell(possibleErrors);
        if(myrank == 0) 
        {
	    cout << "Possible 1 time kmer file size for Rank 0: " << peFileSize / (1024*1024) << "MB " << peFileSize / (Kmer::numBytes()+4) << " kmers." << endl;
        }

        fseek(possibleErrors, 0, SEEK_SET);
        while(ki.read(possibleErrors)) {
                ki.includeCount();
        }

	// discard possible error file (hopefully nothing was actually flushed to disk...)
        unlink(pefilename);
        fseek(possibleErrors, 0, SEEK_SET);
	ftruncate(fileno(possibleErrors), 0);
        fclose(possibleErrors);
        possibleErrors = NULL;

	freeBuffer(peBuffer);
	peBuffer = NULL;

        MPI_Pcontrol(-1,"ReprocessFirstObservations");

        MPI_Pcontrol(1,"HashClean");

        int64_t maxcount = 0;
        int64_t globalmaxcount;
        int64_t hashsize = 0;
        for(auto itr = kmercounts.begin(); itr != kmercounts.end(); ++itr)
        {
		#ifdef KHASH
			if(!itr.isfilled()) continue;   // not all entries are full in khash
			int allcount = get<2>(itr.value());
		#else
			int allcount =  get<2>(itr->second);
		#endif
			if(allcount > maxcount)  maxcount = allcount;
			nonerrorkmers += allcount;
			++hashsize;
        }


        int64_t totalnonerror;
        int64_t distinctnonerror;
        CHECK_MPI( MPI_Reduce(&nonerrorkmers, &totalnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
        CHECK_MPI( MPI_Reduce(&hashsize, &distinctnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
        CHECK_MPI( MPI_Allreduce(&maxcount, &globalmaxcount, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD) );
        if(myrank == 0) 
        {
			cout << "PASS "<< pass << " finished " << endl;
			cout << "Kmerscount hash includes " << distinctnonerror << " distinct elements" << endl;
			cout << "Kmerscount non error kmers count is " << totalnonerror << endl;
			cout << "Global max count is " << globalmaxcount << endl;
			cout << "Large count histogram is of size " << HIGH_NUM_BINS << endl;
        }
#ifdef HISTOGRAM
        vector<int64_t> hist(COUNT_THRESHOLD,0); 	// zero initialize
        vector<int64_t> hist_high(HIGH_NUM_BINS,0);
#endif
        
#ifdef HEAVYHITTERS
	double shh = MPI_Wtime();
        MPI_Op ufxreducempiop;
        CHECK_MPI( MPI_Op_create(MPI_UFXReduce, true, &ufxreducempiop) );    // create commutative mpi-reducer
        CHECK_MPI( MPI_Allreduce(MPI_IN_PLACE, Frequents, nfreqs, MPIType<UFX2ReduceObj>(), ufxreducempiop, MPI_COMM_WORLD) );
        int heavyitems = 0;
        int millioncaught = 0;
        int64_t heavycounts = 0;
        for(size_t i=0; i<nfreqs; ++i)
        {
			Kmer kmer = heavyhitters.Get(i);
			uint64_t myhash = kmer.hash();
			double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
			size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());
			if(owner == myrank)
			{
				MERARR karr = kmer.getArray();
				if(Frequents[i].count)
					kmercounts.insert(make_pair(karr, make_tuple(Frequents[i].ACGTleft,Frequents[i].ACGTrigh,Frequents[i].count)));
				else
					if(myrank  == 0) cout << "Zero count frequent item found, not inserting" << endl;
				nonerrorkmers += Frequents[i].count;
				heavycounts += Frequents[i].count;
				if(Frequents[i].count > maxcount)  maxcount = Frequents[i].count;
				if(Frequents[i].count > MILLION)  ++millioncaught;
				++hashsize;
				++heavyitems;
			}
        }
        delete [] Frequents;    // free memory
        heavyhitters.Clear();
        int totalheavyitems, totalmillioncaught;
        int64_t totalheavycounts;
        CHECK_MPI( MPI_Reduce(&millioncaught, &totalmillioncaught, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) );
        CHECK_MPI( MPI_Reduce(&heavyitems, &totalheavyitems, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) );
        CHECK_MPI( MPI_Reduce(&heavycounts, &totalheavycounts, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
        CHECK_MPI( MPI_Reduce(&nonerrorkmers, &totalnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
        CHECK_MPI( MPI_Reduce(&hashsize, &distinctnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
        CHECK_MPI( MPI_Allreduce(&maxcount, &globalmaxcount, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD) );
	double ehh = MPI_Wtime();
        if(myrank == 0)
        {
			cout << "****** After high freq k-mers are communicated (" << (ehh-shh) << " s)*******" << endl;
			cout << "Total heavy hitters " << totalheavyitems << " for a total of " << totalheavycounts << " counts" << endl;
			cout << "Averaging " << static_cast<double>(totalheavycounts) / static_cast<double>(totalheavyitems) << " count per item" << endl;
			cout << "which caught " << totalmillioncaught << " entries of over million counts" << endl;
			cout << "Kmerscount hash includes " << distinctnonerror << " distinct elements" << endl;
			cout << "Kmerscount non error kmers count is " << totalnonerror << endl;
			cout << "Global max count is " << globalmaxcount << endl;
        }
#endif

        nonerrorkmers = 0;	// reset
        int64_t overonecount = 0;
        auto itr = kmercounts.begin();
        while(itr != kmercounts.end())
        {
#ifdef KHASH
		if(!itr.isfilled()) { ++itr; continue; }
		int allcount = get<2>(itr.value());
#else
		int allcount =  get<2>(itr->second);
#endif
#ifdef HISTOGRAM
		if(allcount <= 0)
			++(hist[0]);
		else if(allcount <= COUNT_THRESHOLD)
			++(hist[allcount-1]);
		else if(allcount <= COUNT_THRESHOLD_HIGH)
			++(hist_high[(allcount - COUNT_THRESHOLD-1)/HIGH_BIN]);
		else
			++(hist_high[HIGH_NUM_BINS-1]);
#endif
		if(allcount < ERR_THRESHOLD)
		{
			--hashsize;
			auto newitr = itr;
			++newitr;
			kmercounts.erase(itr);	// amortized constant	
			// Iterators, pointers and references referring to elements removed by the function are invalidated.
			// All other iterators, pointers and references keep their validity.
			itr = newitr;
		}
		else	
		{
			nonerrorkmers += allcount;
			++itr;
		}
		if(allcount > 1)
		{
			overonecount += allcount;
		}
        }
#ifdef HISTOGRAM
	hist[0] = kmersprocessed - overonecount;	// over-ride the value for 1, which wasn't correct anyway
#endif

        CHECK_MPI( MPI_Reduce(&nonerrorkmers, &totalnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) ); 
        CHECK_MPI( MPI_Reduce(&hashsize, &distinctnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) ); 
        if(myrank == 0) 
        {
			cout << "Erroneous count < " << ERR_THRESHOLD  << " cases removed " << endl;
			cout << "Kmerscount hash includes " << distinctnonerror << " distinct elements" << endl;
			cout << "Kmerscount non error kmers count is " << totalnonerror << endl;
	}
#ifdef HISTOGRAM
	assert( hist.size() == COUNT_THRESHOLD );
	assert( ((int64_t*) hist.data()) + COUNT_THRESHOLD - 1 == &(hist.back()) );
	assert( hist_high.size() == HIGH_NUM_BINS );
	assert( ((int64_t*) hist_high.data()) + HIGH_NUM_BINS - 1 == &(hist_high.back()) );
        if(myrank == 0) 
	{
			CHECK_MPI( MPI_Reduce(MPI_IN_PLACE, hist.data(), COUNT_THRESHOLD, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );
			CHECK_MPI( MPI_Reduce(MPI_IN_PLACE, hist_high.data(), HIGH_NUM_BINS, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );

			stringstream ss, ss1,ss2;
			ss << KMER_LENGTH;
			string hname = "histogram_k";
			hname += ss.str();
			string hnamehigh = hname + "_beyond";
			ss1 << COUNT_THRESHOLD;
			hnamehigh += ss1.str();
			hnamehigh += "_binsize";
			ss2 << HIGH_BIN;
			hnamehigh += ss2.str();
			hname += ".txt";
			hnamehigh += ".txt";
			hname = getRankPath(hname.c_str(), -1);
			hnamehigh = getRankPath(hnamehigh.c_str(), -1);
			ofstream hout(hname.c_str());
			int bin = 0;
			for(auto it = hist.begin(); it != hist.end(); it++) {
				hout << ++bin << "\t" << *it << "\n";
			}
			
			ofstream hhigh(hnamehigh.c_str());
			bin = COUNT_THRESHOLD;
			for(auto it = hist_high.begin(); it != hist_high.end(); it++) {
				hout << bin << "\t" << *it << "\n";
				bin += HIGH_BIN;
			}
        }
        else
        {
			CHECK_MPI( MPI_Reduce(hist.data(), NULL, COUNT_THRESHOLD, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );	// receive buffer not significant at non-root
			CHECK_MPI( MPI_Reduce(hist_high.data(), NULL, HIGH_NUM_BINS, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD) );
        }
	if (myrank == 0) cout << "Generated histograms" << endl;
#endif
        MPI_Pcontrol(-1,"HashClean");
    } // end of if(pass==1)
}


int main(int argc, char ** argv)
{
    CHECK_MPI( MPI_Init(&argc, &argv) );
    CHECK_MPI( MPI_Comm_size(MPI_COMM_WORLD,&nprocs) );
    CHECK_MPI( MPI_Comm_rank(MPI_COMM_WORLD,&myrank) );
	
    double start_time = MPI_Wtime();
	
    nonerrorkmers = 0;
    totaltime = 0;
    kmersprocessed = 0;
    readsprocessed = 0;

    char option;
    bool opt_err = false;
    char *input_fofn = NULL;
    bool use_shm = false;
    string prefix = "";
    bool use_hll = false;

	while ((option = getopt(argc, argv, "Q:f:e:p:sH")) != -1) {
		switch (option) {
		case 'Q':
			phred_encoding = atoi(optarg);
			ignored_ext = phred_encoding + 1;// too low to be a good extension
			break;
		case 'f':
			input_fofn = strdup(optarg);
			break;
		case 'e':
			ERR_THRESHOLD = atoi(optarg);
			break;
		case 'p':
			prefix = optarg;
			break;
		case 's':
			use_shm = true;
			break;
                case 'H':
                        use_hll = true;
                        break;
		default:
			opt_err = true;
		}
	}

	if (opt_err || !input_fofn)
    {
        if(myrank  == 0)
        {
            cout << "Usage: ./ufx -Q phred_encoding -f filename -e error_threshold <-p prefix> <-s> <-H>" << endl;
			cout << "'phred_encoding' is the offset for the quality scores in the fastq file (33 or 64)" << endl;
            cout << "'filename' is a text file containing the paths to each file to be counted - one on each line -" << endl;
            cout << "'error_threshold' is the lowest number of occurrences of a k-mer for it to be considered a real (non-error) k-mer" << endl;
            cout << "'prefix' is optional and if set, the code prints all k-mers with the prefix" << endl;
            cout << "-s is optional and if set, shared memory will be used to cache all files" << endl;
            cout << "-H is optional and if set, then HyperLogLog will be used to estimate the required bloomfilter" << endl;
        }
        MPI_Finalize();
        return 0;
    }

	if(myrank  == 0)
	{
		cout << "HipMER k-mer characterization (ufx generation) step" << endl;
		cout << "You are running with the following settings:" << endl;
		cout << "K-mer length = " << KMER_LENGTH << endl;
		cout << "Max k-mer length (internal) = " << MAX_KMER_SIZE << endl;
		cout << "Phred encoding = " << phred_encoding << endl;
	}

	double t01 = MPI_Wtime();

	vector<filedata> allfiles = GetFiles(input_fofn);
	Kmer::set_k(KMER_LENGTH);
	double cardinality;
		//use_hll标志表示 是否使用HyperLogLog算法来估计KMER数
        if (use_hll) {
   		ProudlyParallelCardinalityEstimate(allfiles, cardinality, use_shm);	// doesn't update kmersprocessed yet (but updates readsprocessed)
        } else {
		// just estimate the cardinality using fastq size and sampling.
		long mySums[3] = {0,0,0};
		long &myCardinality = mySums[0];
		long &myTotalReads = mySums[1];
		long &myTotalBases = mySums[2];
		for(int i = 0; i < allfiles.size(); i++) {
			if ((i+1) % nprocs == myrank) {
				long totalReads, totalBases;
				long fileSize = estimate_fq(allfiles[i].filename, 5000, &totalReads, &totalBases);
				myTotalReads += totalReads;
				myTotalBases += totalBases;
				assert(fileSize == allfiles[i].filesize);
				if (totalReads > 0) {
					long kmersPerRead = ((totalBases+totalReads-1) / totalReads) - KMER_LENGTH + 1;
					myCardinality += kmersPerRead * totalReads;
					cout << "Cardinality for " << allfiles[i].filename << ": " << myCardinality << endl;
				}
			}
		}
		MPI_Allreduce(MPI_IN_PLACE, mySums, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
		if (myrank == 0) {
			cout << "Estimated our cardinality: " << myCardinality << " totalReads: " << myTotalReads << " totalBases: " << myTotalBases << endl;
		}
		// baseline of 10M kmers
		if (myCardinality < 10000000) myCardinality = 10000000;
		// assume maximum of 90% of the kmers are unique, because at least some have to repeat...
		cardinality = 0.9 * myCardinality / (double) nprocs;
	}
#ifndef HIPMER_BLOOM64
	if (cardinality > 1L<<31) {
		cout << "Reduced cardinality to fit within int32_t" << endl;
		cardinality = 1L<<31 - 1L;
	}
#endif
	double t02 = MPI_Wtime();
	if (myrank == 0) {
       		cout << "Stage 1 - Cardinality: " << t02-t01 << endl;
	}

        // only one pass is necessary now
	ProcessFiles(allfiles, 1, cardinality, use_shm);	// determine final hash-table entries using bloom filter
	double t03 = MPI_Wtime();
	if (myrank == 0) {
       		cout << "Stage 2 - Exchange and Count: " << t03-t02 << endl;
	}

	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
	double t04 = MPI_Wtime();
	if(myrank  == 0)
	{
       		cout << "Stage 3 (noop): " << t04-t03 << endl;
	}
	
	int64_t sendbuf = kmercounts.size(); 
	int64_t recvbuf, totcount, maxcount;
	CHECK_MPI( MPI_Exscan(&sendbuf, &recvbuf, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Allreduce(&sendbuf, &totcount, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD) );
	CHECK_MPI( MPI_Allreduce(&sendbuf, &maxcount, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD) );
	int64_t totkmersprocessed;
	CHECK_MPI( MPI_Reduce(&kmersprocessed, &totkmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD) );
	if(myrank  == 0)
	{
		double imbalance = static_cast<double>(nprocs * maxcount) / static_cast<double>(totcount);  
		cout << "Load imbalance for final k-mer counts is " << imbalance << endl;
		cout << "CardinalityEstimate " << static_cast<double>(totkmersprocessed) / (MEGA * (t02-t01) * nprocs) << " MEGA k-mers per sec/proc in " << (t02-t01) << " seconds"<< endl;
		cout << "Bloom filter + hash pre-allocation  " << static_cast<double>(totkmersprocessed) / (MEGA * (t03-t02) * nprocs) << " MEGA k-mers per sec/proc in " << (t03-t02) << " seconds" << endl;
		cout << "Actual hash-based counting  " << static_cast<double>(totkmersprocessed) / (MEGA * (t04-t03) * nprocs) << " MEGA k-mers per sec/proc in " << (t04-t03) << " seconds" << endl;
	}
	CHECK_MPI( MPI_Barrier(MPI_COMM_WORLD) );
	double t05 = MPI_Wtime();
	if (myrank == 0)
	{
		cout << "Stage 4 - loadImbalance: " << t05-t04 << endl;
	}
	vector<filedata>().swap(allfiles);
	string countfile = string(input_fofn)+".ufx.bin";


#ifndef BENCHMARKONLY
	if(myrank  == 0) cout << "Writing to binary via MPI-IO" << endl;
    
	MPI_File thefile;
	int64_t lengthuntil;

	if (!use_shm) {
		CHECK_MPI( MPI_File_open(MPI_COMM_WORLD, (char*) countfile.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile) );
		if(myrank == 0)	lengthuntil = 0;
		else	lengthuntil = recvbuf;
	}

#ifdef DEBUG
    
    char mergraphfilepath[MAX_FILE_PATH];
    sprintf(mergraphfilepath, "%s.mergraph.%d", countfile.c_str(), myrank);
    char *mergraphfile = get_rank_path(mergraphfilepath, myrank);
    ofstream mergraph(mergraphfile);
#endif

	int ufxcount = kmercounts.size();
	ufxpack * packed = new ufxpack[ufxcount];

	int packedIdx = 0;
	for(auto itr = kmercounts.begin(); itr != kmercounts.end(); ++itr)
	{
    #ifdef KHASH
		if(!itr.isfilled()) { continue; }
		auto key = itr.key();
		auto val = itr.value();
    #else
		auto key = itr->first;
		auto val = itr->second;;
    #endif

		Kmer kmer(key);
		ufxpack upack;
		upack.arr = key;
		PackIntoUFX(get<0>(val), get<1>(val), get<2>(val), upack);
		packed[packedIdx++] = upack;
#ifdef DEBUG
        mergraph << kmer.toString() << "\t" << get<0>(val)[0] << "\t" << get<0>(val)[1]<< "\t" << get<0>(val)[2]<< "\t" << get<0>(val)[3] << "\t0\t" << (get<2>(val)-get<0>(val)[0]-get<0>(val)[1]-get<0>(val)[2]-get<0>(val)[3]) << "\t" << get<1>(val)[0] << "\t" << get<1>(val)[1] << "\t" << get<1>(val)[2] << "\t" << get<1>(val)[3] << "\t0\t" << (get<2>(val)-get<1>(val)[0]-get<1>(val)[1]-get<1>(val)[2]-get<1>(val)[3]) << "\t0" << endl;
#endif
	}
#ifndef KHASH   // otherwise the destructor will probably get to it
	kmercounts.clear();
#endif
    
	if (!use_shm) {
		MPI_Datatype datatype;
		CHECK_MPI( MPI_Type_contiguous(sizeof(ufxpack), MPI_CHAR, &datatype ) );
		CHECK_MPI( MPI_Type_commit(&datatype) );
		int dsize;
		CHECK_MPI( MPI_Type_size(datatype, &dsize) );

		int mpi_err = MPI_File_set_view(thefile, lengthuntil * dsize, datatype, datatype, (char*)"external32", MPI_INFO_NULL);
		if (mpi_err == 51) {
			// external32 datarep is not supported, use native instead
			CHECK_MPI( MPI_File_set_view(thefile, lengthuntil * dsize, datatype, datatype, (char*)"native", MPI_INFO_NULL) );
		} else {
			CHECK_MPI(mpi_err);
		}

		MPI_Status status;
		CHECK_MPI( MPI_File_write_all(thefile, packed, ufxcount, datatype, &status) );
		CHECK_MPI( MPI_File_close(&thefile) );
	} else {
		stringstream ss;
		string rank;
		ss << myrank;
		ss >> rank;
		unsigned found = countfile.find_last_of("/");
		countfile = countfile.substr(found+1);
		countfile += rank;
		string sizefile = countfile + string(".entries");	// filename<pid>.entries
		int64_t ufxcount64bit = static_cast<int64_t>(ufxcount);
		int64_t bytes2write =  ufxcount64bit * static_cast<int64_t>(sizeof(ufxpack));
    
		int fd = shm_open(countfile.c_str(), O_CREAT | O_TRUNC | O_RDWR, S_IRUSR | S_IWUSR); // O_TRUNC: if file exists, the file is truncated to zero size
		if (fd <= 0) {
			printf("Error, shm_open failed %d %s when trying to open %s\n", errno, strerror(errno), countfile.c_str());
			exit(-1);
		}
		ftruncate(fd, bytes2write);
		void * tmpaddr = mmap(NULL, bytes2write , PROT_READ | PROT_WRITE , MAP_SHARED , fd , 0);
		if (tmpaddr == MAP_FAILED) {
			cout << "mmap open for write failed with " << errno << " - " << strerror(errno) << endl;
		}
		else {
			if(myrank == 0) cout << "Opened mmap file " << countfile << " ufx binary" << endl;
		}
		memcpy(tmpaddr, packed, bytes2write);
		close(fd);

		/* Now write the file that has the total size */
		int fs = shm_open(sizefile.c_str(), O_CREAT | O_TRUNC | O_RDWR, S_IRUSR | S_IWUSR); // O_TRUNC: if file exists, the file is truncated to zero size
		if (fs <= 0) {
			printf("Error, shm_open failed %d %s when trying to open %s\n", errno, strerror(errno), sizefile.c_str());
			exit(-1);
		}
		ftruncate(fs, sizeof(int64_t));
		int64_t totalufxcount;
		CHECK_MPI( MPI_Allreduce(&ufxcount64bit, &totalufxcount, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD) );
		void * tmp2add = mmap(NULL, sizeof(int64_t) , PROT_READ | PROT_WRITE , MAP_SHARED , fs , 0);	// now write the total size
		if (tmp2add == MAP_FAILED) {
			cout << "mmap open for write failed with " << errno << " - " << strerror(errno) << endl;
		}
		else {
			if(myrank == 0) cout << "Opened mmap file " << sizefile << " that lists total ufx entries " << totalufxcount << endl;
		}
		memcpy(tmp2add, &totalufxcount, sizeof(int64_t));
		close(fs);
    
	}

	double t06 = MPI_Wtime();
	if (myrank == 0)
	{
		cout << "Stage 6 - FileIO: " << t06-t05 << endl;
	}

	if(myrank == 0)	cout << "File write completed\n";
	if (!use_shm) 
	{
		if(prefix != "")
		{
			string sanity = "ufx_";
			stringstream ss;
			ss << myrank;
			sanity += ss.str();
			sanity += string(".");
			sanity += prefix;
			ofstream outsanity(sanity.c_str());

			int prelen = prefix.length();
			for(int i=0; i < ufxcount; ++i)
			{
				Kmer kmer(packed[i].arr);
				string kstr = kmer.toString();
				if(kstr.substr(0, prelen) == prefix)
				{
					outsanity << kstr << "\t" << packed[i].count << "\t" << packed[i].left << "[" << packed[i].leftmin << "," << packed[i].leftmax;
					outsanity << "]" << packed[i].right << "[" << packed[i].rightmin << "," << packed[i].rightmax << "]" << endl;
				}
			}
			outsanity.close();
		}
		
		if(myrank  == 0)
		{
			cout << "Finished writing, here is the top of processor 0's data" << endl;
			for(int i=0; i< 10 && i < ufxcount; ++i)
			{
				Kmer kmer(packed[i].arr);
				cout << kmer.toString() << "\t" << packed[i].count << "\t" << packed[i].left << "[" << packed[i].leftmin << "," << packed[i].leftmax;
				cout << "]" << packed[i].right << "[" << packed[i].rightmin << "," << packed[i].rightmax << "]" << endl;
			}
			cout << "Here is a random sample from the global data" << endl;
#ifdef DEBUG__DISABLED
#if GCC_VERSION < 40300
			auto generator = tr1::mt19937();
			tr1::uniform_int<int64_t> distribution(0,count-1);
			generator.seed((unsigned int)time(NULL));
#else
			auto generator = mt19937();
			uniform_int_distribution<int64_t> distribution(0,count-1);
			generator.seed((unsigned int)time(NULL));
#endif
				
			FILE * f = fopen(countfile.c_str(), "r");
			if(!f)
			{
				cerr << "Problem reading binary input file\n";
				return 1;
			}	
			for(int i=0; i< 10; ++i)
			{
				int64_t randindex  = distribution(generator);
				ufxpack upack;
				fseek (f, randindex * static_cast<int64_t>(dsize), SEEK_SET );
				fread(&upack, dsize, 1, f);
				Kmer kmer(upack.arr);
				cout << "LOC " << randindex << ":\t" << kmer.toString() << "\t\t";
				cout << upack.left << "[" << upack.leftmin << "," << upack.leftmax << "]";
				cout << upack.right << "[" << upack.rightmin << "," << upack.rightmax << "]" << endl;
			} 
#endif // #ifdef DEBUG__DISABLED
		}
	} // use_shm
  	delete [] packed;
#endif // #ifndef BENCHMARKONLY
	double t07 = MPI_Wtime();
	if (myrank==0) 
	{
		cout << "Stage 7 - shmem storage: " << t07-t06 << endl;
	}

	double elapsed_time = MPI_Wtime() - start_time;
	double slowest_time;
	CHECK_MPI( MPI_Allreduce(&elapsed_time, &slowest_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) );
	if (!myrank)
		printf("Overall time for %s is %.2f s\n", basename(argv[0]), slowest_time);

   	MPI_Finalize();
    return 0;
}

