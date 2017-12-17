#!/bin/bash
# Set SGE options:
#$ -clear
#$ -l exclusive.c
#$ -cwd
#$ -j y
#$ -pe pe_16 16
#$ -w e

set -e

N=${N:=16}
if [ -n "${NSLOTS}" ]
then
  N=${NSLOTS}
fi

if [ -n "$1" ]
then
  HIPMER_INSTALL=$1
else
  HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-genepool}
fi

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
  
  export CORES_PER_NODE=16
  module list
  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-genepool-validation-${THREADS}-${JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir -p $RUNDIR
  echo "Copying validate set to ${RUNDIR}"
  rsync -av ${INST}/etc/meraculous/pipeline/ ${RUNDIR}/

  echo "Running with ${THREADS} threads!"
  
  MPIRUN="mpirun -n" UPCRUN="upcrun -n" DEBUG=1 UPC_PTHREADS=0 ${INST}/bin/run_hipmer_validation.sh ${INST} ${RUNDIR}/meraculous-validation.config 

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

