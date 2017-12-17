#!/bin/bash
#PBS -q debug
#PBS -l mppwidth=24
#PBS -j oe

N=$(cat ${PBS_NODEFILE} | wc -l)

HIPMER_INSTALL=${HIPMER_INSTALL:=${SCRATCH}/install-hopper}
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
  export RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-hopper-ecoli-diploid-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
  mkdir $RUNDIR

  HIPMER_ECOLI_DATA=${HIPMER_ECOLI_DATA:=${SCRATCH}/hipmer_ecoli_data}
  HIPMER_ECOLI_DATA_URL=http://portal.nersc.gov/archive/home/r/regan/www/HipMer/hipmer_ecoli.tar.gz
  if [ ! -d ${HIPMER_ECOLI_DATA} ]
  then
    echo "Downloading ecoli fastq from the NERSC HPSS system... this could take a while"
    if [ ${N} -gt 0 ]
    then
      echo "Do not download the data within a batch job. run 'mkdir ${HIPMER_ECOLI_DATA}; lfs setstripe -c 72 ${HIPMER_ECOLI_DATA} ; cd ${HIPMER_ECOLI_DATA} ; curl --max-time 7200 -o - $HIPMER_ECOLI_DATA_URL | tar -xzf' offline"
      exit 1
    fi
    (mkdir ${HIPMER_ECOLI_DATA}; lfs setstripe -c 72 ${HIPMER_ECOLI_DATA} ; cd ${HIPMER_ECOLI_DATA} ; curl --max-time 7200 -o - $HIPMER_ECOLI_DATA_URL | tar -xzf -)
    [ -f ${HIPMER_ECOLI_DATA}/ECO.fastq.0.fq ] && [ -f ${HIPMER_ECOLI_DATA}/ECO.fastq.1.fq ]
  fi
  for i in ${HIPMER_ECOLI_DATA}/* ; do ln -s $i $RUNDIR ; done
  export ILLUMINA_VERSION=13

  cp -p ${INST}//etc/meraculous/pipeline/meraculous-ecoli-diploid.config $RUNDIR 

  UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB:=1500}
  MPIRUN="aprun -n" UPCRUN="upcrun -n" UPC_SHARED_HEAP_MB=${UPC_SHARED_HEAP_MB} ${INST}/bin/run_hipmer.sh $INST ${RUNDIR}/meraculous-ecoli-diploid.config

else
  echo "Could not find HipMer installed at ${INST}"
  exit 1
fi

