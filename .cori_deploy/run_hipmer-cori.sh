#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00
#SBATCH --job-name=HipMer

set -x
set -e


# evaluate any key=value pairs on the command line
for arg in "$@"
do
  if [ "${arg/=}" != "${arg}" ]
  then
    eval "export ${arg%=*}='${arg##*=}'"
    echo "placed into the env: ${arg}"
    shift
  else
    break
  fi
done

export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}
N=${N:=${SLURM_NTASKS}}

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-cori}
if [ -f "$1" ]
then
  HIMPER_CONFIG=$1
  shift
fi
if [ -d "$1" ]
then
  HIPMER_INSTALL=${1}
  shift
fi

HIPMER_CONFIG=$(readlink -f ${HIPMER_CONFIG})
USAGE="$0 [ ENV=val ... ] [ meraculous.config ] [ hipmer-install-path ]
or set environment:
  HIPMER_INSTALL=${HIPMER_INSTALL}
  HIPMER_CONFIG=${HIPMER_CONFIG}
  RUNDIR=${RUNDIR}
"

if [ ! -f "${HIPMER_CONFIG}" ]
then
  echo "$USAGE"
  echo "No config found: $HIPMER_CONFIG"
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
  
  module list
  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-${HIPMER_CONFIG##*/}-${THREADS}-${SLURM_JOB_NAME}-${SLURM_JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir -p $RUNDIR

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=3000}
  cd $RUNDIR
  MPIRUN="srun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${HIPMER_INSTALL}/bin/run_hipmer.sh ${HIPMER_INSTALL} ${HIPMER_CONFIG}

else
  echo "Could not find HipMer installed at ${HIPMER_INSTALL}"
  exit 1
fi

