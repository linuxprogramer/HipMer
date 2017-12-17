#!/bin/bash

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-generic}
USAGE="$0 hipmer-install-path
or set environment: HIPMER_INSTALL=${HIPMER_INSTALL}
"

if [ -z "$CORES_PER_NODE" ]; then
	export CORES_PER_NODE=`grep processor /proc/cpuinfo|wc -l`
fi

N=${CORES_PER_NODE}
export CAN_SPLIT_JOB=1

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
  
  export THREADS=${THREADS:=${N}}
  #export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-generic-validation-${THREADS}-$(date '+%Y%m%d_%H%M%S')}
  export RUNDIR=${RUNDIR:=$SCRATCH/validation-${THREADS}-$(date '+%Y%m%d_%H%M%S')}
  mkdir $RUNDIR

  HIPMER_VALIDATION_DATA=${HIPMER_VALIDATION_DATA:=${SCRATCH}/hipmer_validation_data}
  if [ ! -d ${HIPMER_VALIDATION_DATA} ]
  then
      mkdir ${HIPMER_VALIDATION_DATA}
      echo "Copying validate set to ${HIPMER_VALIDATION_DATA}"
      rsync -av ${INST}/etc/meraculous/pipeline/ ${HIPMER_VALIDATION_DATA}/
      [ -f ${HIPMER_VALIDATION_DATA}/frags.fastq.25K ] && [ -f ${HIPMER_VALIDATION_DATA}/jumps.fastq.25K ]
  fi
  for i in ${HIPMER_VALIDATION_DATA}/* ; do ln -s $i $RUNDIR ; done

  cp -p ${INST}//etc/meraculous/pipeline/meraculous-validation.config $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="mpirun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} POST_RUN="check_validation_results.sh" \
    ${INST}/bin/run_hipmer.sh ${INST} ${RUNDIR}/meraculous-validation.config 

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

