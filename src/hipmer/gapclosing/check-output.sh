#!/bin/bash

ddir=$1
reference_final=$2

function get_field {
    echo -n "$1   "
    grep "$1" $ddir/merauder_*.log|$HOME/code/meraculous/gapclosing/avg.py $2 2
}    

function get_timings {
    x=`grep Elapsed $ddir/merauder_*.log|grep $1|$HOME/code/meraculous/gapclosing/avg.py $2 2`
	if [ "$x" != "" ]; then
		echo "$1  $x"
	fi
}    

if [ -f "$ddir/assembly.0.fa" ]; then
    # get_timings scaffold 7
    # get_timings "  srf-counts" 5
    # get_timings "  srf-reads" 5
    # get_timings "  calc" 7
    # get_timings contigs 6
    # get_timings blastmap 5
    # get_timings "  direct" 7
    # get_timings "  projected" 7
    # get_timings "  get_alignment" 5
    # get_timings "  put_read_orient" 5
    # get_timings "  remainder" 5
    # get_timings printing 5
    # get_timings FASTQ 6
    # get_timings "  reading" 6
    # get_timings "    get_read_orient" 5
    # get_timings "    gap_copy" 5
    # get_timings "  atomics" 5
    # get_timings "closing" 5
    # get_timings "  splinting" 6
    # get_timings "  walks" 7
    # get_timings "  patching" 6
    # get_timings "FASTA" 4
    # get_timings "thread " 5

    # get_field "Peak memory" 3

    if [ $reference_final ]; then
        echo -n "differences "
        cat $ddir/assembly.*.fa |sort > lala; diff lala $reference_final |wc -l
    fi
fi
