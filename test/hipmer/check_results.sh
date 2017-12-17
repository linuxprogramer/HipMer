#!/bin/bash

if [ "$RUNDIR" != "" ]; then
    pushd $RUNDIR > /dev/null
fi
INUFX=ALL_INPUTS.fofn.ufx.bin
if [ -z "$THREADS" ]; then
	export THREADS=`grep processor /proc/cpuinfo|wc -l`
fi

MPIRUN=${MPIRUN:="mpirun -n $THREADS"}
if [ ! -f "$INUFX.txt.full.sorted" ]; then
    [ -f "${INUFX}.txt.full" ] || $MPIRUN b2tufx-${KMER_LENGTH} ${INUFX} >> b2tufx.log 2>&1
    if [ ! -f "${INUFX}.txt.full" ]; then
      echo  "No ${INUFX}.txt.full file!" 1>&2
    else 

      nlines=$(wc -l ${INUFX}.txt.full|cut -f1 -d ' ')
      if (( nlines < 100 )); then
        cat ${INUFX}.txt.full
        exit 1
      fi
      sort $INUFX.txt.full > $INUFX.txt.full.sorted 
    fi
fi

[ -f canonical_contigs.fa ] || cat output_*_contigs.fasta | canonical_assembly.pl > canonical_contigs.fa

PASS="\x1B[92mPASS\x1B[0m"
FAIL="\x1B[91mFAIL\x1B[0m"
echo "Evaluating checksums:"
echo "$ufx_md5sum  $INUFX.txt.full.sorted" |md5sum --check &> /dev/null
if [ "$?" == 0 ]; then
    echo -e "  UFX:                " $PASS
else
    echo -e "  UFX:                " $FAIL
fi
echo "$contigs_md5sum  canonical_contigs.fa" |md5sum --check &> /dev/null
if [ "$?" == 0 ]; then
    echo -e "  contigs:            " $PASS
else
    echo -e "  contigs:            " $FAIL
fi
echo "$assembly_md5sum  final_assembly.fa" |md5sum --check &> /dev/null
if [ "$?" == 0 ]; then
    echo -e "  assembly:           " $PASS
else
    echo -e "  assembly:           " $FAIL
    tr '[:lower:]' '[:upper:]' < final_assembly.fa > final_assembly.uppercase.fa
    echo "$uc_assem_md5sum  final_assembly.uppercase.fa" |md5sum --check &> /dev/null
    if [ "$?" == 0 ]; then
        echo -e "  uppercase assembly: " $PASS
    else
        echo -e "  uppercase assembly: " $FAIL
    fi
fi
if [ "$RUNDIR" != "" ]; then
    popd > /dev/null
fi

