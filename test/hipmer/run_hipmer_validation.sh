#!/bin/bash
set -e

# This script copies the validation datasets and configs and executes the hipmer pipline
USAGE="$0 install-path meraculous.config"
INSTALL_PATH=$(readlink -f $1)

if [ ! -d "${INSTALL_PATH}" ] || [ ! -f $INSTALL_PATH/etc/meraculous/pipeline/frags.fastq.25K ]
then
  echo "$USAGE" 1>&2
  exit 1
fi

RUNDIR=${RUNDIR:=$(hostname)-$$-$(date '+%Y%m%d-%H%M%S')}
mkdir -p $RUNDIR
export RUNDIR=$RUNDIR
cp -p $INSTALL_PATH/etc/meraculous/pipeline/frags.fastq.25K $INSTALL_PATH/etc/meraculous/pipeline/jumps.fastq.25K $RUNDIR/

ret=0
export DEBUG=1 # make sure ufx.bin.txt is created
export PATH=${INSTALL_PATH}/bin:${PATH}
run_hipmer.sh $@ || ret=1

cd $RUNDIR
if [ ! -f ALL_INPUTS.fofn.ufx.bin ]
then
  echo "Could not find ALL_INPUTS.fofn.ufx.bin!"
  exit 1
fi

if [ -f ALL_INPUTS.fofn.ufx.bin.txt.full ]
then
  [ -f ALL_INPUTS.fofn.ufx.bin.txt.full.sorted ] || sort ALL_INPUTS.fofn.ufx.bin.txt.full > ALL_INPUTS.fofn.ufx.bin.txt.full.sorted 
fi

[ -f canonical_contigs.fa ] || cat output_*_contigs.fasta | canonical_assembly.pl  > canonical_contigs.fa

check_validation_results.sh || ret=$?

exit $ret

