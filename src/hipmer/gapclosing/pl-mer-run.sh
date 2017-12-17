#!/bin/bash

genome=$1
gaps=$2

if [ "$genome" == "ecoli" ]; then
    datadir="../data/ecoli"
    ./merauder4.pl -i216 -m31 -D3 -R2.75 -Q64 -A -P \
        -s $datadir/p3.1.srf -c $datadir/contigs.fa -g $datadir/$gaps 
elif [ "$genome" == "human" ]; then
    datadir="../data/human/"
    ./merauder4.pl -i407 -m51 -D4 -P -R1.75 -Q33 -A -V \
        -s $datadir/p3.srf -c $datadir/hs.m51.bubbletigs.fa -g $datadir/$gaps 
else 
    echo "Invalid genome: ecoli or human"
fi





