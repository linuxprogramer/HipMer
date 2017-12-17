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
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "Deleter.h"
#include "Friends.h"

using namespace std;
template <class T>
bool from_string(T & t, const string& s, std::ios_base& (*f)(std::ios_base&))
{
   istringstream iss(s);
   return !(iss >> f >> t).fail();
}

std::string itos(int n)
{
   const int max_size = std::numeric_limits<int>::digits10 + 1 /*sign*/ + 1 /*0-terminator*/;
   char buffer[max_size] = {0};
   sprintf(buffer, "%d", n);
   return std::string(buffer);
}
		
int myrank;

int main(int argc, char ** argv)
{	
	int nprocs; 
    	MPI_Init(&argc, &argv);
    	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	
	MPI_Datatype datatype;
        MPI_Type_contiguous(sizeof(kmerpack), MPI_CHAR, &datatype );
        MPI_Type_commit(&datatype);
        int dsize;
        MPI_Type_size(datatype, &dsize);

	int64_t entries2dump;

	struct stat filestatus;
  	stat( argv[ 1 ], &filestatus );
  	cout << filestatus.st_size << " bytes\n";
	int64_t fileentries = static_cast<int64_t>(filestatus.st_size) / static_cast<int64_t>(dsize);
	int64_t threshold = 0;
	if(argc < 3)
	{
		entries2dump = fileentries;

	}
	else
	{
		int64_t n;
		from_string(n,string(argv[2]),std::dec);

		if(fileentries > n)
		{
			entries2dump = n;
		}
		else
		{
			entries2dump = fileentries;	// all we have
		}

		if(argc > 3)
		{
			from_string(threshold,string(argv[3]),std::dec);
			cout << "Dumping only entries above " << threshold;
		}
	}
	Kmer::set_k(KMER_LENGTH);

	int64_t perproc = entries2dump / static_cast<int64_t>(nprocs);
	int64_t perme;
	if(myrank == nprocs-1)	perme = entries2dump - perproc*(nprocs-1);
	else	perme = perproc;

  	cout << "dumping " << perme << " entries\n";

	int64_t mybegin = perproc * static_cast<int64_t>(myrank) ;
	FILE * f = fopen(argv[1], "r");
	if(!f)
	{
		cerr << "Problem reading binary input file\n";
		return 1;
	}
		
	kmerpack kpack;
	fseek (f, mybegin * static_cast<int64_t>(dsize), SEEK_SET );

	string outname = string(argv[1])+ itos(myrank)+string(".txt");
	FILE * pFile = fopen (outname.c_str(),"w");
	for(int64_t i=0; i< perme; ++i)
	{
		fread(&kpack, dsize, 1, f);
		Kmer kmer(kpack.arr);
		if(kpack.count > threshold)
			fprintf(pFile, "%s\t\t%d\n", kmer.toString().c_str(), kpack.count);
	}
	fclose(f);
	fclose(pFile);
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank == 0)
	{
		string command =  "cat";
        	for(int i=0; i<nprocs; ++i)
        	{
                	command += " ";
			string name = string(argv[1])+ itos(i)+string(".txt");
                	command += name;
		}
        	command += " > ";
		string finalname = string(argv[1])+ string(".txt.full");
        	command += finalname;
        	cout << command << endl;
        	system(command.c_str());
	}

	MPI_Finalize();
	return 0;
}
