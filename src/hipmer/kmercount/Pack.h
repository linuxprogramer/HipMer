#ifndef _PACK_H_
#define _PACK_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
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
#include <locale>
#include "HipMERdefines.h"
#include "Kmer.hpp"
#include "Friends.h"
#include "FriendsMPI.h"

using namespace std;

extern int nprocs;
extern int myrank;


#ifdef HEAVYHITTERS
#ifndef MAXHITTERS
#define MAXHITTERS 200000
#endif
extern SimpleCount<Kmer> heavyhitters;
extern  UFX2ReduceObj * Frequents;
#endif


// The bloom filter pass; extensions are ignored
inline size_t FinishPackPass1(vector< vector<Kmer> > & outgoing, Kmer & kmerreal)
{
    uint64_t myhash = kmerreal.hash();  // whichever one is the representative
    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());
    
#ifdef HEAVYHITTERS
    if(heavyhitters.IsMember(kmerreal))
    {
        // no-op here
        // cout << kmerreal.toString() << " is high freq" << endl;
    }
    else {
#endif
        outgoing[owner].push_back(kmerreal);
#ifdef HEAVYHITTERS
    }
#endif
    return outgoing[owner].size();
}

inline void IncrementBasedOnQual(array<int,4> & mergraphentry, char qual, char seq)
{
    if(qual >= phred_encoding+extensionCutoff)
    {
        switch(seq)
        {
            case 'A' :
                mergraphentry[0]++;
                break;
            case 'C' :
                mergraphentry[1]++;
                break;
            case 'G' :
                mergraphentry[2]++;
                break;
            case 'T' :
                mergraphentry[3]++;
                break;
        }
    }
}

// The hash table pass; extensions are important
inline size_t FinishPackPass2(vector< vector<Kmer> > & outgoing, vector<vector<array<char,2>>> & extquals,
                      vector<vector<array<char,2>>> & extseqs, Kmer & kmerreal, array<char,2> & extqual, array<char,2> & extseq)
{
    uint64_t myhash = kmerreal.hash();  // whichever one is the representative
    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());

#ifdef HEAVYHITTERS
    size_t location = heavyhitters.FindIndex(kmerreal);    // maxsie if not found
    if(location < heavyhitters.maxsize)
    {
        Frequents[location].count += 1; // add one more count
        IncrementBasedOnQual(Frequents[location].ACGTleft, extqual[0], extseq[0]);
        IncrementBasedOnQual(Frequents[location].ACGTrigh, extqual[1], extseq[1]);
    }
    else {
        //cout << myrank << ": " << kmerreal.toString() << " is NOT a heavy hitter " << endl;
#endif
    outgoing[owner].push_back(kmerreal);
    extquals[owner].push_back(extqual);
    extseqs[owner].push_back(extseq);
#ifdef HEAVYHITTERS
    }
#endif
    return outgoing[owner].size();
}


size_t PackEnds(string & seq, string & qual, int j, vector< vector<Kmer> > & outgoing, vector<vector<array<char,2>>> & extquals,
                  vector<vector<array<char,2>>> & extseqs, int pass, int lastCountedBase)
{
    bool isCounted = lastCountedBase >= j + KMER_LENGTH;
    if (pass == 1 && !isCounted) return 0;

    string kmerextstr;
    try {
        kmerextstr = seq.substr(j, KMER_LENGTH);
    } catch (std::exception const &exc) {
        std::cerr << "Exception caught in file " << __FILE__<< " at " << __LINE__ << ": " << exc.what() << "\n";
        std::cerr << "j = " << j << ", KMER_LENGTH = " << KMER_LENGTH << ", len seq = " << seq.length() << "\n";
        std::cerr << "seq = " << seq << "\n";
        MPI_Finalize();
        exit(1);
    }
    for (auto & c: kmerextstr) c = toupper(c);	// convert all to uppercase
    size_t found=kmerextstr.find('N');
    if (found!=std::string::npos) return 0;	// if there is an 'N', toss it

    size_t procSendCount;
    Kmer kmerreal(kmerextstr.c_str());
    if(pass == 1)
    {
        kmerreal = kmerreal.rep();
        procSendCount = FinishPackPass1(outgoing, kmerreal);
    }
    else if(pass == 2)   // otherwise we don't care about the extensions
    {
        bool isLeft = j > 0;
        bool isRight = j < seq.length()-KMER_LENGTH;
        array<char,2> extseq  = {isLeft ? seq[j-1] : '-', isRight ? seq[j+KMER_LENGTH] : '-'};
        array<char,2> extqual = {(char) (isLeft ? (char) qual[j-1] : (char) ignored_ext), (char) (isRight ? (char) qual[j+KMER_LENGTH] : (char) ignored_ext)};
	if (!isCounted) {
		// Do not count this kmer, but extensions may be counted...
		extqual[0] = 0 - extqual[0];
		extqual[1] = 0 - extqual[1];
	}
        Kmer kmertwin = kmerreal.twin();
        if(kmertwin < kmerreal)	// the real k-mer is not lexicographically smaller
        {
            kmerreal = kmertwin;
            RevCompSwitch(extseq[0]);
            RevCompSwitch(extseq[1]);
            std::swap(extseq[0], extseq[1]);	// swap the sequence ends
            std::swap(extqual[0], extqual[1]);	// swap the quality ends
        }
        procSendCount = FinishPackPass2(outgoing, extquals, extseqs, kmerreal, extqual, extseq);
    }
    return procSendCount;
}

#endif
