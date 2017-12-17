#!/bin/bash
#PBS -q debug
#PBS -l mppwidth=24
#PBS -j oe

N=$(cat ${PBS_NODEFILE} | wc -l)

export CORES_PER_NODE=24

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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-chr14-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
  echo "Preparing ${RUNDIR}"
  mkdir -p $RUNDIR
  # to improve IO performance on hopper set stripe to max on scratch
  lfs setstripe -c -1 ${RUNDIR}
  echo "Setting stripe to maximum for $RUNDIR"

  $INST/bin/hipmer_setup_chr14_data.sh ${INST} ${RUNDIR}

#  export ILLUMINA_VERSION=18

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="aprun -n" UPCRUN="upcrun -q -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-chr14.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

