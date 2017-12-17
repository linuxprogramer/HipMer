#!/bin/bash

USAGE="$0 HIPMER_INSTALL RUNDIR
 (or set them as environmental variables)
 HIPMER_INSTALL=/path/to/hipmer
 RUNDIR=/path/to/desired/working_directory
 SCRATCH=/path/to/rundir (if RUNDIR is not specified)
 THREADS=concurrency (if not set by PBS_NODEFILE automatically)

"
set -x
set -e

INST=${HIPMER_INSTALL:=$1}
SCRATCH=${SCRATCH:=/tmp}
if [ -n "$2" ]
then
  RUNDIR=${RUNDIR:=$2}
else
  RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-$(hostname)-ecoli-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
fi

export HIPMER_INSTALL=${INST}
export RUNDIR

if [ -d ${INST}/bin ] && [ -x ${INST}/bin/run_hipmer.sh ] && [ -n "${RUNDIR}" ]
then
  echo "Setting up ${RUNDIR} with E. coli data set at $(date)"
else
  echo "$USAGE"
  exit 1
fi

mkdir -p $RUNDIR

echo "Linking E. coli set to ${RUNDIR}"
cp -p ${INST}/etc/meraculous/pipeline/*ecoli*.config ${RUNDIR}/


HIPMER_ECOLI_DATA=${HIPMER_ECOLI_DATA:=${SCRATCH}/hipmer_ecoli_data}
HIPMER_ECOLI_DATA_URL=http://portal.nersc.gov/archive/home/r/regan/www/HipMer/hipmer_ecoli.tar.gz
if [ ! -d ${HIPMER_ECOLI_DATA} ]
then
    echo "Downloading ecoli fastq from the NERSC HPSS system... this could take a while"
    if [  -n "$PBS_NODEFILE" ]
    then
      echo "Do not download the data within a batch job. run 'mkdir ${HIPMER_ECOLI_DATA}; lfs setstripe -c 72 ${HIPMER_ECOLI_DATA} ; cd ${HIPMER_ECOLI_DATA} ; curl --max-time 7200 -o - $HIPMER_ECOLI_DATA_URL | tar -xzf' offline"
      exit 1
    fi
    (mkdir ${HIPMER_ECOLI_DATA}; lfs setstripe -c 72 ${HIPMER_ECOLI_DATA} || /bin/true ; cd ${HIPMER_ECOLI_DATA} ; curl --max-time 7200 -o - $HIPMER_ECOLI_DATA_URL | tar -xzf -)
    [ -f ${HIPMER_ECOLI_DATA}/ECO.fastq.0.fq ] && [ -f ${HIPMER_ECOLI_DATA}/ECO.fastq.1.fq ]
fi
for i in ${HIPMER_ECOLI_DATA}/* ; do ln -f $i $RUNDIR || ln -fs $i $RUNDIR ; done

cp -p ${INST}/etc/meraculous/pipeline/meraculous-ecoli*.config $RUNDIR 


