#!/bin/bash
# Set SGE options:
#$ -cwd
#$ -j y
#$ -pe pe_16 16
#$ -l ram.c=7.5G

set -x
set -e

N=16
if [ -n "${NSLOTS}" ]
then
  N=${NSLOTS}
fi


HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-genepool}
if [ -f "$1" ]
then
  HIPMER_CONFIG=${1}
  shift
fi
if [ -d "${1}" ]
then
  HIPMER_INSTALL=${1}
fi

USAGE="$0 [meraculous.config] [HIPMER_INSTALL]
or set environment:
  HIPMER_INSTALL=${HIPMER_INSTALL}
  HIPMER_CONFIG=${HIPMER_CONFIG}
  RUNDIR=${RUNDIR}
"

if [ ! -f "${HIPMER_CONFIG}" ]
then
  echo "$USAGE"
  exit 1
fi

if [ -d ${HIPMER_INSTALL}/bin ] && [ -x ${HIPMER_INSTALL}/bin/run_hipmer.sh ]
then
  src=$(echo ${HIPMER_INSTALL}/env*.sh)
  if [ -f ${src} ]
  then
    . ${src}
  else
    echo "Could not find an environment file to source in ${HIPMER_INSTALL}!" 1>&2
  fi
  
  module list || /bin/true
  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$(pwd)/hipmer-genepool-${HIPMER_CONFIG##*/}-${THREADS}-${JOB_ID}-$(date '+%Y%m%d_%H%M%S')}

  echo "Preparing ${RUNDIR}"
  mkdir -p $RUNDIR

  echo "Running with ${THREADS} threads!"
  
  MPIRUN="mpirun -n" UPCRUN="upcrun -n" ${HIPMER_INSTALL}/bin/run_hipmer.sh $HIPMER_INSTALL ${HIPMER_CONFIG}

else
  echo "Could not find HipMer installed at ${HIPMER_INSTALL}"
  exit 1
fi


