#!/bin/bash

genome=$1
threads=$2

root_ddir="."


if [ "$genome" == "ecoli" ]; then
    datadir="$root_ddir/ecoli"
    rm -rf $datadir/results/2015-*
    upcrun -q -shared-heap=1G -n$threads\
        ./merauder -i216 -I10 -m19 -D4 -R2.75 -Q64 -A \
        -f $datadir/ECO.fastq.\*.fq \
		-s $datadir/FsrfFile_ \
        -c $datadir/FASTAcontigs_ \
        -a $datadir/merAlignerOutput_ \
		-o $datadir/results \
		-b 0,1
    rm -rf $datadir/results/latest
    ln -s `ls $datadir/results/|grep 2015` $datadir/results/latest
    ./check-output.sh $datadir/results/latest $datadir/results/assembly-all-sorted
elif [ "$genome" == "human" ]; then
    datadir="$root_ddir/human"
	upcrun -q -shared-heap=1GB -n$threads \
		./merauder.$threads -i393 -I43 -m51 -D4 -R1.75 -A -P \
		-f $datadir/s_1_\*_sequence.fq \
        -s $datadir$threads/FsrfFile_ \
        -c $datadir$threads/FASTAcontigs_ \
        -a $datadir$threads/merAlignerOutput_ \
		-o $datadir/results/n$threads \
        -b 0,1
elif [ "$genome" == "human-chr14" ]; then
    datadir="$root_ddir/human-chr14"
    rm -rf $datadir/results/2015-*
	upcrun -q -shared-heap=2GB -n$threads \
		./merauder -i159 -I17 -m51 -D4 -R1.75 -A -P \
		-f $datadir/frag_\*fastq \
        -s $datadir/n$threads/FsrfFile_ \
        -c $datadir/n$threads/FASTAcontigs_ \
        -a /dev/shm/merAlignerOutput_ \
		-o $datadir/results/ \
        -b 0,1
        #-a $datadir/n$threads/merAlignerOutput_ 
    rm -rf $datadir/results/latest
    ln -s `ls $datadir/results/|grep 2015` $datadir/results/latest
    ./check-output.sh $datadir/results/latest $datadir/results/assembly-all-sorted
elif [ "$genome" == "wheat" ]; then
    datadir="$root_ddir/wheat"
	upcrun -q -shared-heap=1500MB -n$threads \
		./merauder.$threads -i736 -I26 -m51 -D3 -R1.75 -A \
		-f $datadir/7\*.fastq.lane\* \
		-s $datadir/run$threads/RsrfFile_ \
        -c $datadir/FASTAcontigs_ \
        -a $datadir/merAlignerOutput_ \
		-o $datadir/results/n$threads \
	    -b 0,1,2,3,4,5,10,11,12,13,14,15
else 
    echo "Invalid genome: ecoli, human or wheat"
fi






