#!/bin/bash

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-generic}
USAGE="$0 hipmer-install-path
or set environment: HIPMER_INSTALL=${HIPMER_INSTALL}
"
if [ -z "$CORES_PER_NODE" ]; then
	export CORES_PER_NODE=`grep processor /proc/cpuinfo|wc -l`
fi

N=${CORES_PER_NODE}
export CAN_SPLIT_JOB=1

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
  
  export THREADS=${THREADS:=${N}}

  export RUNDIR=${RUNDIR:=$SCRATCH/chr14-${THREADS}-$(date '+%Y%m%d_%H%M%S')}
  mkdir $RUNDIR

  HIPMER_CHR14_DATA=${HIPMER_CHR14_DATA:=${SCRATCH}/hipmer_chr14_data}
  HIPMER_CHR14_DATA_URL=http://portal.nersc.gov/archive/home/s/shofmeyr/www/hipmer/chr14.tar.gz
  if [ ! -d ${HIPMER_CHR14_DATA} ]
  then
    echo "Downloading chr14 fastq... this could take a while"
    (mkdir ${HIPMER_CHR14_DATA}; cd ${HIPMER_CHR14_DATA} ; curl --max-time 7200 -o - $HIPMER_CHR14_DATA_URL | tar -xzf -)
    [ -f ${HIPMER_CHR14_DATA}/frag_1.fastq ] && [ -f ${HIPMER_CHR14_DATA}/frag_2.fastq ]
  fi
  for i in ${HIPMER_CHR14_DATA}/* ; do ln -s $i $RUNDIR ; done
#  export ILLUMINA_VERSION=18

  cp -p ${INST}//etc/meraculous/pipeline/meraculous-chr14.config $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="mpirun -n" UPCRUN="upcrun -q -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-chr14.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

