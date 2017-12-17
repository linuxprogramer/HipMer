#!/bin/bash
#PBS -q debug
#PBS -l mppwidth=24
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
  
  module list
  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-$(hostname)-validation-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
  echo "Preparing ${RUNDIR}"
  mkdir -p $RUNDIR

  ${INST}/bin/hipmer_setup_validation_data.sh ${INST} ${RUNDIR}

  echo "Running with ${THREADS} threads!"
  
  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="aprun -n" UPCRUN="upcrun -n" DEBUG=1 UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} POST_RUN=check_validation_results.sh \
    ${INST}/bin/run_hipmer_validation.sh ${INST} ${RUNDIR}/meraculous-validation.config 

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

