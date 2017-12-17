#!/bin/bash
# Set SGE options:
#$ -clear
#$ -l exclusive.c
#$ -cwd
#$ -j y
#$ -pe pe_16 64

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
  module list || /bin/true
  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-genepool-chr14-${THREADS}-${JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir -p $RUNDIR

  echo "Preparing ${RUNDIR}"
  mkdir -p $RUNDIR

  $INST/bin/hipmer_setup_chr14_data.sh ${INST} ${RUNDIR}

  echo "Running with ${THREADS} threads!"
  
  MPIRUN="mpirun -n" UPCRUN="upcrun -n" ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-chr14.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi


