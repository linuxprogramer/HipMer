#!/bin/bash

# Compares results for two runs on the same genome. Because of non-determinancy in the code, the
# final assemblies will likely differ, so this test compares a number of statistics for equivalence.

dir1=$1
dir2=$2

pushd $dir1 > /dev/null
lib_name1=`ls *.fofn|grep -v ALL_INPUTS|cut -f1 -d'.'`
popd > /dev/null
pushd $dir2 > /dev/null
lib_name2=`ls *.fofn|grep -v ALL_INPUTS|cut -f1 -d'.'`
popd > /dev/null
if [ "$lib_name1" != "$lib_name2" ]; then
	echo "DIFFER: Libraries are different $lib_name1 != $lib_name2"
fi

# won't work for shm
# diff -qs $dir1/ALL_INPUTS.fofn.ufx.bin $dir2/ALL_INPUTS.fofn.ufx.bin

identical_tests=0
tot_tests=0

for f in bmaMeta-splinter "$lib_name1-bmaMeta-spans" "$lib_name2-insert.txt" \
    nContigs_bubble_contigs.fasta.txt nUUtigs_contigs.fasta.txt nGaps_SsrfFile-1.txt \
    nLinks-1.txt nScaffolds_SsrfFile-1.txt nUUtigs_contigs.fasta.txt Round1-bmaMeta; do
	if [ -f "$dir1/$f" ]; then
		if [ ! -f "$dir2/$f" ]; then
			echo DIFFER: $f found in $dir1 but not $dir2
		else
			x=`diff -qs $dir1/$f $dir2/$f|grep identical`
			if [ "$x" != "" ]; then
				echo "[$tot_tests] Identical: $f"
				identical_tests=$((identical_tests+1))
			else
				echo "[$tot_tests] DIFFER: $f"
				diff -W40 -y $dir1/$f $dir2/$f
			fi
		fi
	    tot_tests=$((tot_tests+1))
    fi
done

assembly1=`wc -l $dir1/final_assembly.fa|cut -f1 -d' '`
assembly2=`wc -l $dir2/final_assembly.fa|cut -f1 -d' '`

if [ "$assembly1" != "$assembly2" ]; then
	echo "[$tot_tests] DIFFER: counts in final_assembly.fa: $assembly1 != $assembly2"
else
	echo "[$tot_tests] Identical: counts in final_assembly.fa"
	identical_tests=$((identical_tests+1))
fi
tot_tests=$((tot_tests+1))

gstats1=`grep "Main genome" $dir1/gapclosing.log`
gstats2=`grep "Main genome" $dir2/gapclosing.log`

if [ "$gstats1" != "$gstats2" ]; then
	echo "[$tot_tests] DIFFER: stats for final assembly:"
	echo "  -----------"
	echo "$gstats1"
	echo "  -----------"
	echo "$gstats2"
	echo "  -----------"
else
	echo "[$tot_tests] Identical: stats for final_assembly.fa"
	identical_tests=$((identical_tests+1))
fi
tot_tests=$((tot_tests+1))

echo "-> Identical in $identical_tests/$tot_tests tests"

