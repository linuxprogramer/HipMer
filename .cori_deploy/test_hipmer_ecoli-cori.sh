#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00
#SBATCH --job-name=HipmerEcoli

set -x
set -e

export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}
N=${N:=${SLURM_NTASKS}}

if [ -d "$1" ]
then
  HIPMER_INSTALL=${1}
fi
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
  
  module list
  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$SCRATCH/ecoli-${THREADS}-${SLURM_JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir -p $RUNDIR

  $INST/bin/hipmer_setup_ecoli_data.sh ${INST} ${RUNDIR}

  cp -p ${INST}/etc/meraculous/pipeline/meraculous-ecoli.config $RUNDIR 
  
  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=3000}
  MPIRUN="srun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} POST_RUN=check_ecoli_results.sh \
    ${INST}/bin/run_hipmer.sh ${INST} ${RUNDIR}/meraculous-ecoli.config 

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

