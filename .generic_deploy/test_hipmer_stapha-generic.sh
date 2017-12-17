#!/bin/bash

#export USE_SHM=1

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

  export RUNDIR=${RUNDIR:=$SCRATCH/stapha-${THREADS}-$(date '+%Y%m%d_%H%M%S')}
  mkdir $RUNDIR

  HIPMER_STAPHA_DATA=${HIPMER_STAPHA_DATA:=${SCRATCH}/hipmer_stapha_data}
  HIPMER_STAPHA_DATA_URL=http://portal.nersc.gov/archive/home/s/shofmeyr/www/hipmer/stapha.tar.gz
  if [ ! -d ${HIPMER_STAPHA_DATA} ]
  then
    echo "Downloading stapha fastq... this could take a while"
    (mkdir ${HIPMER_STAPHA_DATA}; cd ${HIPMER_STAPHA_DATA} ; curl --max-time 7200 -o - $HIPMER_STAPHA_DATA_URL | tar -xzf -)
    [ -f ${HIPMER_STAPHA_DATA}/frag_1.fastq ] && [ -f ${HIPMER_STAPHA_DATA}/frag_2.fastq ]
  fi
  for i in ${HIPMER_STAPHA_DATA}/* ; do ln -s $i $RUNDIR ; done
#  export ILLUMINA_VERSION=18

  cp -p ${INST}//etc/meraculous/pipeline/meraculous-stapha.config $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="mpirun -n" UPCRUN="upcrun -q -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-stapha.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

