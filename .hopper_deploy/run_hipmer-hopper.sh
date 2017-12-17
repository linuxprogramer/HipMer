#!/bin/bash
#PBS -q debug
#PBS -l mppwidth=24
#PBS -j oe

cd $PBS_O_WORKDIR
N=$(cat ${PBS_NODEFILE} | wc -l)

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-hopper}
CONFIG=${CONFIG:=$(echo *.config)}
USAGE="$0 hipmer-install-path
or set environment:
  HIPMER_INSTALL=${HIPMER_INSTALL}
  CONFIG=${CONFIG}
  ILLUMINA_VERSION=${ILLUMINA_VERSION}
  UPC_SHARED_HEAP_MB
"
if [ ! -f ${CONFIG} ]
then
  echo "Please set the CONFIG environmental variable"
  exit 1
fi
CONFIG=$(readlink -f ${CONFIG})

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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-${CONFIG##*/}-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir $RUNDIR

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=900}
  cd $RUNDIR
  MPIRUN="aprun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${CONFIG}

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

