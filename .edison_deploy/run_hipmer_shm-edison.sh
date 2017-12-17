#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --job-name=HipMer

export USE_SHM=1
export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}
N=${N:=${SLURM_NTASKS}}

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-edison}
USAGE="$0
set environment:
 HIPMER_INSTALL=${HIPMER_INSTALL}
 HIPMER_CONFIG=${HIPMER_CONFIG}
"
INST=${HIPMER_INSTALL}
if [ ! -f "${HIPMER_CONFIG}" ]
then
  echo "Please specify a valid HIPMER_CONFIG in the environment"
  exit 1
fi

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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-shm-${HIPMER_CONFIG##*/}-${THREADS}-${SLURM_JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  echo "Preparing ${RUNDIR}"
  mkdir -p $RUNDIR

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="srun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh ${INST} ${HIPMER_CONFIG}

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

