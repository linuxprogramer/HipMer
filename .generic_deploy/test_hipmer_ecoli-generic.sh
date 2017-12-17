#!/bin/bash

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-generic}
USAGE="$0 hipmer-install-path
or set environment: HIPMER_INSTALL=${HIPMER_INSTALL}
"

CONFIG_FILE=meraculous-ecoli.config
#CONFIG_FILE=meraculous-ecoli-diploid.config

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

  export RUNDIR=${RUNDIR:=$SCRATCH/ecoli-${THREADS}-$(date '+%Y%m%d_%H%M%S')}
  mkdir $RUNDIR

  HIPMER_ECOLI_DATA=${HIPMER_ECOLI_DATA:=${SCRATCH}/hipmer_ecoli_data}
  HIPMER_ECOLI_DATA_URL=http://portal.nersc.gov/archive/home/r/regan/www/HipMer/hipmer_ecoli.tar.gz
  if [ ! -d ${HIPMER_ECOLI_DATA} ]
  then
    echo "Downloading ecoli fastq from the NERSC HPSS system... this could take a while"
    (mkdir ${HIPMER_ECOLI_DATA}; cd ${HIPMER_ECOLI_DATA} ; curl --max-time 7200 -o - $HIPMER_ECOLI_DATA_URL | tar -xzf -)
    [ -f ${HIPMER_ECOLI_DATA}/ECO.fastq.0.fq ] && [ -f ${HIPMER_ECOLI_DATA}/ECO.fastq.1.fq ]
  fi
  for i in ${HIPMER_ECOLI_DATA}/* ; do ln -s $i $RUNDIR ; done

  cp -p ${INST}//etc/meraculous/pipeline/$CONFIG_FILE $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="mpirun -n" UPCRUN="upcrun -q -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} POST_RUN="check_ecoli_results.sh" \
    ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/$CONFIG_FILE

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

