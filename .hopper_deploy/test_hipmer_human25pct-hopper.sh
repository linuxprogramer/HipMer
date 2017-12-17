#!/bin/bash
#PBS -q debug
#PBS -l mppwidth=600
#PBS -l walltime=30:00
#PBS -j oe

set -x
set -e

N=$(cat ${PBS_NODEFILE} | wc -l)

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-hopper}
USAGE="$0 hipmer-install-path
or set environment: HIPMER_INSTALL=${HIPMER_INSTALL}
"
INST=${HIPMER_INSTALL:=$1}
if [ -d ${INST}/bin ] && [ -x ${INST}/bin/run_hipmer.sh ]
then
  src=$(echo ${INST}/env*.sh)
  if [ -f ${src} ]
  then
    . ${src}
  else
    echo "Could not find an environment file to source in ${INST}!" 1>&2
  fi
  
  module list || /bin/true

  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-hopper-human25pct-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir $RUNDIR

  HIPMER_HUMAN_DATA=${HIPMER_HUMAN_DATA:=${SCRATCH}/hipmer_human_data}
  HIPMER_HUMAN_DATA_URL=http://portal.nersc.gov/archive/home/r/regan/www/HipMer/hipmer_human_subsets.tar.gz
  if [ ! -d ${HIPMER_HUMAN_DATA} ] && [ -f ${HIPMER_HUMAN_DATA}/bwa-25pct-read1.fastq ]
  then
    echo "Downloading ecoli fastq from the NERSC HPSS system... this could take a while"
    if [ ${N} -gt 0 ]
    then
      echo "Do not download this data in a batch job.  Do it offline:  mkdir ${HIPMER_HUMAN_DATA}; lfs setstripe -c 72 ${HIPMER_HUMAN_DATA} ; cd ${HIPMER_HUMAN_DATA} ; curl --max-time 7200 -o - $HIPMER_HUMAN_DATA_URL | tar -xzf -"
      exit 1
    fi
    (mkdir -p ${HIPMER_HUMAN_DATA}; lfs setstripe -c 72 ${HIPMER_HUMAN_DATA} ; cd ${HIPMER_HUMAN_DATA} ; curl --max-time 7200 -o - $HIPMER_HUMAN_DATA_URL | tar -xzf -)
    [ -f ${HIPMER_HUMAN_DATA}/bwa-25pct-read1.fastq ] && [ -f ${HIPMER_HUMAN_DATA}/bwa-25pct-read2.fastq ]
  fi
  for i in ${HIPMER_HUMAN_DATA}/*25pct*.fastq ; do ln -s $i $RUNDIR ; done
  export ILLUMINA_VERSION=18

  cp -p ${INST}/etc/meraculous/pipeline/meraculous-human25pct.config $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="aprun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-human25pct.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi


