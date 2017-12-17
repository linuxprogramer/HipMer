#!/bin/bash

datadir="/work/shofmeyr/meraculous/e2/"

./gapPlacer.pl -m 31 -i 216:10 \
-b $datadir/combined-clean.bm \
-s $datadir/combined.srf \
-f $datadir/ECO.fastq.info.0 \
-c $datadir/combined-clean.fa


