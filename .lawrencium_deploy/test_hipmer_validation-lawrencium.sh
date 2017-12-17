#!/bin/bash -l
#SBATCH --job-name=HipmerValidation
#SBATCH --partition=lr3
#SBATCH --qos=lr_debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=0:10:0

set -x
set -e

export GASNET_BACKTRACE=1
export OMPI_MCA_mpi_warn_on_fork=0
#export CAN_SPLIT_JOB=0
env

export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}
N=${N:=${SLURM_NTASKS}}
echo "Detected CORES_PER_NODE=${CORES_PER_NODE} and N=${N}"

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-lawrencium}
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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-validation-${THREADS}-${SLURM_JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  echo "Preparing ${RUNDIR}"
  mkdir -p ${RUNDIR}
  ${INST}/bin/hipmer_setup_validation_data.sh ${INST} ${RUNDIR}
  
  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="mpirun -n" UPCRUN="upcrun -n" DEBUG=1 UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} POST_RUN=check_validation_results.sh \
     ${INST}/bin/run_hipmer_validation.sh ${INST} ${RUNDIR}/meraculous-validation.config 

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

