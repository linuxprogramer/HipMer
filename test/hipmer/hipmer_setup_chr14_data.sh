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
  RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-$(hostname)-chr14-${THREADS}-${PBS_JOBID}-$(date '+%Y%m%d_%H%M%S')}
fi

export HIPMER_INSTALL=${INST}
export RUNDIR

if [ -d ${INST}/bin ] && [ -x ${INST}/bin/run_hipmer.sh ] && [ -n "${RUNDIR}" ]
then
  echo "Setting up ${RUNDIR} with human chr14 data set at $(date)"
else
  echo "$USAGE"
  exit 1
fi

HIPMER_CHR14_DATA=${HIPMER_CHR14_DATA:=${SCRATCH}/hipmer_chr14_data}
HIPMER_CHR14_DATA_URL=http://portal.nersc.gov/archive/home/s/shofmeyr/www/hipmer/chr14.tar.gz

if [ ! -d ${HIPMER_CHR14_DATA} ]
then
    echo "Downloading human chr14 fastq... this could take a while"
    if [ -n "${PBS_NODEFILE}" ]
    then
      echo "Do not download this data in a batch job.  Do it offline:  mkdir ${HIPMER_CHR14_DATA}; lfs setstripe -c 72 ${HIPMER_CHR14_DATA} ; cd ${HIPMER_CHR14_DATA} ; curl --max-time 7200 -o - $HIPMER_CHR14_DATA_URL | tar -xzf -"
      exit 1
    fi
    (mkdir ${HIPMER_CHR14_DATA}; lfs setstripe -c 72 ${HIPMER_CHR14_DATA} || /bin/true ; cd ${HIPMER_CHR14_DATA} ; curl --max-time 7200 -o - $HIPMER_CHR14_DATA_URL | tar -xzf -)
    [ -f ${HIPMER_CHR14_DATA}/frag_1.fastq ] && [ -f ${HIPMER_CHR14_DATA}/frag_2.fastq ]
fi
for i in ${HIPMER_CHR14_DATA}/* ; do ln -f $i $RUNDIR || ln -fs $i $RUNDIR ; done

cp -p ${INST}/etc/meraculous/pipeline/meraculous-chr14.config $RUNDIR 



