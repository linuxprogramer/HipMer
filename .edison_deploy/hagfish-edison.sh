#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=480
#PBS -l walltime=1:00:00
#PBS -j oe

set -x
set -e

N=$(cat ${PBS_NODEFILE} | wc -l)

min_cores=480

if [ "$N" -lt "$min_cores" ]; then
	echo "Too few cores, $N, for this dataset. Rerun with at least $min_cores cores."
	exit 1
fi

export CORES_PER_NODE=24

HIPMER_HAGFISH_DATA=/scratch3/scratchdirs/shofmeyr/hagfish/HAGFISH

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-edison}
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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-hagfish-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
  echo "Preparing ${RUNDIR}"
  mkdir $RUNDIR
  # to improve IO performance on edison set stripe to max on scratch
  lfs setstripe -c -1 ${RUNDIR}
  echo "Setting stripe to maximum for $RUNDIR"

  for i in ${HIPMER_HAGFISH_DATA}/* ; do ln -s $i $RUNDIR ; done

  cp -p ${INST}/etc/meraculous/pipeline/meraculous-hagfish.config $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="aprun -n" UPCRUN="upcrun -q -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-hagfish.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi


