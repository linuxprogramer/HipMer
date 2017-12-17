#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00
#SBATCH --job-name=HipMer

export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}
N=${N:=${SLURM_NTASKS}}

set -x
set -e


HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-cori}
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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-human25pct-${THREADS}-${SLURM_JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  echo "Preparing ${RUNDIR}"
  mkdir -p $RUNDIR

  $INST/bin/hipmer_setup_human_subsets_data.sh ${INST} ${RUNDIR}

  export ILLUMINA_VERSION=18

  cp -p ${INST}/etc/meraculous/pipeline/meraculous-human25pct.config $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="srun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-human25pct.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi


