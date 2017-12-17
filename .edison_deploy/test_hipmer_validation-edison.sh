#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:30:00
#SBATCH --job-name=HipMer

export CORES_PER_NODE=${CORES_PER_NODE:=${SLURM_TASKS_PER_NODE%%\(*}}
N=${N:=${SLURM_NTASKS}}

set -x
set -e

export CORES_PER_NODE=24

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
  
  module list
  export THREADS=${THREADS:=${N}}
  export RUNDIR=${RUNDIR:=$SCRATCH/validation-${THREADS}-${SLURM_JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir -p $RUNDIR

  HIPMER_VALIDATION_DATA=${HIPMER_VALIDATION_DATA:=${SCRATCH}/hipmer_validation_data}
  if [ ! -d ${HIPMER_VALIDATION_DATA} ]
  then
      mkdir ${HIPMER_VALIDATION_DATA}
      echo "Copying validate set to ${HIPMER_VALIDATION_DATA}"
      rsync -av ${INST}/etc/meraculous/pipeline/ ${HIPMER_VALIDATION_DATA}/
      [ -f ${HIPMER_VALIDATION_DATA}/frags.fastq.25K ] && [ -f ${HIPMER_VALIDATION_DATA}/jumps.fastq.25K ]
  fi
  for i in ${HIPMER_VALIDATION_DATA}/* ; do [ -f $RUNDIR/${i##*/} ] || ln -s $i $RUNDIR ; done

  cp -p ${INST}/etc/meraculous/pipeline/meraculous-validation.config $RUNDIR 
  
  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="srun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} POST_RUN=check_validation_results.sh \
    ${INST}/bin/run_hipmer.sh ${INST} ${RUNDIR}/meraculous-validation.config 

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

