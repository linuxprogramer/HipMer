#!/bin/bash
# Set SGE options:
#$ -clear
#$ -l exclusive.c
#$ -cwd
#$ -j y
#$ -pe pe_16 16

set -x
set -e

N=16
if [ -n "${NSLOTS}" ]
then
  N=${NSLOTS}
fi

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-genepool}
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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-genepool-ecoli-${THREADS}-${JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir -p $RUNDIR
  echo "Copying ecoli set to ${RUNDIR}"
  $INST/bin/hipmer_setup_ecoli_data.sh ${INST} ${RUNDIR}

  echo "Running with ${THREADS} threads!"
  
  MPIRUN="mpirun -n" UPCRUN="upcrun -n" DEBUG=1 UPC_PTHREADS=0 POST_RUN=check_ecoli_results.sh \
    ${INST}/bin/run_hipmer.sh ${INST} ${RUNDIR}/meraculous-ecoli.config 

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

