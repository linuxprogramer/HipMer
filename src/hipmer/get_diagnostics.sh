#!/bin/bash

cd per_thread

echo -e "\n==============\nUFX\n==============\n"
egrep "Kmerscount|Erroneous" ufx.log
echo ""

echo -e "\n==============\ncontigs\n==============\n"
echo -n "Number of UU tigs: "
cat nUUtigs_contigs.fasta.txt
echo ""

echo -e "\n==============\ncontigMerDepth\n==============\n"
grep "Number of records" contigMerDepth.log
echo ""


echo -e "\n==============\ncontigEndAnalyzer\n==============\n"
grep "Number of records" contigEndAnalyzer.log
echo ""

if [ -f bubbleFinder.log ]; then
    echo -e "\n==============\nbubbleFinder\n==============\n"
    egrep "Total bubbles|Total diplot|Total iso" bubbleFinder.log
    echo -n "Number of bubble contigs: " 
    cat nContigs_bubble_contigs.fasta.txt
    echo ""
fi

echo -e "\n==============\ncontigMerDepth2\n==============\n"
grep "Number of records" contigMerDepth2.log
echo ""


echo -e "\n==============\nmerAligner\n==============\n"
for f in `ls merAligner-*.log`; do
    name=${f%.log}; echo ${name#merAligner-}
    grep -A7 "Total reads processed" $f
    echo ""
done


echo -e "\n==============\nhistogrammer\n==============\n"
for f in `ls histogrammer-*.log`; do
    name=${f%.log}; echo ${name#histogrammer-}
    grep "Total blast map lines" $f
    grep -A3 "Hits ommitted" $f
    egrep "Carried out|non-standard|very different" $f
    echo ""
done

echo -e "\n==============\nsplinter\n==============\n"
for f in `ls splinter-*.log`; do
    name=${f%.log}; echo ${name#splinter-}
    grep -A4 "Splints found" $f
    echo ""
done

echo -e "\n==============\nspanner\n==============\n"
for f in `ls spanner-*.log`; do
    name=${f%.log}
    lib=${name#spanner-}
    lib=${lib%-*}
    echo $lib
    echo -n "Total alignments: "
    cat $lib-nTotalAlignments.txt
    grep -A4 "Unused " $f
    grep "Total spans found" $f
    echo ""
done

echo -e "\n==============\nbmaToLinks\n==============\n"
for f in `ls bmaToLinks-*.log`; do 
    grep -A10 "LIB: " $f
    echo -n "Number of links: "
    name=${f%.log}
    round=${name#bmaToLinks-}
    cat nLinks-$round.txt
    echo ""
done


echo -e "\n==============\noNo\n==============\n"
for i in `seq 1 100`; do
    if [ ! -f "oNo-$i-3.log" ]; then
        break
    fi
    bestP=`cat n50_SsrfFile-$i-* |sort -nr -k2 |head -n1|awk '{print $1}'`
    bestN50=`cat n50_SsrfFile-$i-${bestP}.txt|awk '{print $2}'`
    nScaffolds=`cat nScaffolds_SsrfFile-$i-${bestP}.txt`
    onoF=oNo-$i-$bestP.log
    echo "Round $i, selected p $bestP, N50 $bestN50"
    grep "Modal depth" $onoF
    grep -A8 "BEST TIE" $onoF|grep -v "======"
    echo "Number of scaffolds: $nScaffolds"
    echo ""
done

echo -e "\n==============\nGapclosing\n==============\n"
#libs=`grep library gapclosing_0.log|awk -F'[:,]' '{print $2, $3}'`
#for l in libs; do
echo "Unused alignments"
echo -n "    5-TRUNCATED "
grep "5-TRUNC" 0*/gapclosing_*.log|awk '{sum+=$2} END {print sum}'
echo -n "    3-TRUNCATED "
grep "3-TRUNC" 0*/gapclosing_*.log|awk '{sum+=$2} END {print sum}'
echo -n "    NO_ASSOCIATED_GAP "
grep "NO_ASSOC" 0*/gapclosing_*.log|awk '{sum+=$2} END {print sum}'
echo "Total reads placed in gaps"
echo -n "    Aligned "
grep "placed in gaps" 0*/gapclosing_*.log|awk '{sum+=$7} END {print sum}'
echo -n "    Projected "
grep "placed in gaps" 0*/gapclosing_*.log|awk '{sum+=$10} END {print sum}'
echo -n "Reads found "
grep "reads, total" 0*/gapclosing_*.log|awk '{sum+=$5} END {print sum}'
#done
echo -n "Number of spans "
grep Closures 0*/gapclosing_*.log|awk '{sum+=$2} END {print sum}'
echo -n "Number of walks "
grep walks 0*/gapclosing_*.log|awk '{sum+=$4} END {print sum}'
grep "Total gaps closed" gapclosing.log
echo ""

echo -e "\n==============\nStatistics\n==============\n"
echo "Contigs"
echo -en "\t"
grep "Found " upc-canonical-contigs.log
grep "Size" upc-canonical-contigs.log
grep "N50" upc-canonical-contigs.log
echo "Scaffolds"
echo -en "\t"
grep "Found " upc-canonical-assembly.log
grep "Size" upc-canonical-assembly.log
grep "N50" upc-canonical-assembly.log

cd ..