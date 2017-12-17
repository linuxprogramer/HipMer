#!/bin/bash

USAGE="$0 HIPMER_INSTALL RUNDIR
 (or set them as environmental variables)
 HIPMER_INSTALL=/path/to/hipmer
 RUNDIR=/path/to/desired/working_directory
 SCRATCH=/path/to/rundir (if RUNDIR is not specified)
 THREADS=concurrency (if not set by PBS_NODEFILE automatically)

"
set -x
set -e

INST=${HIPMER_INSTALL:=$1}
SCRATCH=${SCRATCH:=/tmp}
if [ -n "$2" ]
then
  RUNDIR=${RUNDIR:=$2}
else
  RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-$(hostname)-validation-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
fi

export HIPMER_INSTALL=${INST}
export RUNDIR

if [ -d ${INST}/bin ] && [ -x ${INST}/bin/run_hipmer.sh ] && [ -n "${RUNDIR}" ]
then
  echo "Setting up ${RUNDIR} with validation data set at $(date)"
else
  echo "$USAGE"
  exit 1
fi

mkdir -p $RUNDIR

echo "Copying validate set to ${RUNDIR}"
cp -p ${INST}/etc/meraculous/pipeline/* ${RUNDIR}/


