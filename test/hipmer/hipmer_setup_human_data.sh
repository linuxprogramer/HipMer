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
  RUNDIR=${RUNDIR:=$SCRATCH/hipmer-test-$(hostname)-human-${THREADS}-${PBS_JOBID}${SLURM_JOB_ID}${JOB_ID}-$(date '+%Y%m%d_%H%M%S')}
fi

export HIPMER_INSTALL=${INST}
export RUNDIR

if [ -d ${INST}/bin ] && [ -x ${INST}/bin/run_hipmer.sh ] && [ -n "${RUNDIR}" ]
then
  echo "Setting up ${RUNDIR} with human data set at $(date)"
else
  echo "$USAGE"
  exit 1
fi

mkdir -p $RUNDIR

HIPMER_HUMAN_DATA=${HIPMER_HUMAN_DATA:=${SCRATCH}/hipmer_human_data}
HIPMER_HUMAN_DATA_URL=http://portal.nersc.gov/archive/home/r/regan/www/HipMer/hipmer_human_fastq.tar.gz 
HIPMER_HUMAN_DATA_URL2=http://portal.nersc.gov/archive/home/r/regan/www/HipMer/hipmer_human_jmp.tar.gz 
HIPMER_HUMAN_DATA_URL3=http://portal.nersc.gov/archive/home/r/regan/www/HipMer/hipmer_human_fos.tar.gz 
if [ ! -d ${HIPMER_HUMAN_DATA} ]
then
    echo "Downloading human fastq from the NERSC HPSS system... this could take a while"
    if [ -n "${PBS_JOBID}" ] || [ -n "${SLURM_JOB_ID}" ] || [ -n "${JOB_ID}" ]
    then
      echo "Do not download this data in a batch job.  Do it offline:  mkdir ${HIPMER_HUMAN_DATA}; lfs setstripe -c 72 ${HIPMER_HUMAN_DATA} ; cd ${HIPMER_HUMAN_DATA} ; curl --max-time 36000 -o - $HIPMER_HUMAN_DATA_URL | tar -xzf -"
      exit 1
    fi
    (mkdir -p ${HIPMER_HUMAN_DATA}.tmp
     lfs setstripe -c 72 ${HIPMER_HUMAN_DATA}.tmp || /bin/true
     cd ${HIPMER_HUMAN_DATA}.tmp && 
       ((curl --max-time 36000 -o - $HIPMER_HUMAN_DATA_URL | tar -xzf -) &
        (curl --max-time 36000 -o - $HIPMER_HUMAN_DATA_URL2 | tar -xzf -) &
        (curl --max-time 36000 -o - $HIPMER_HUMAN_DATA_URL3 | tar -xzf -) &
       wait ) &&
     cd .. && mv ${HIPMER_HUMAN_DATA}.tmp ${HIPMER_HUMAN_DATA}
    )
    [ -f ${HIPMER_HUMAN_DATA}/s_1_1_sequence.permuted ] && [ -f ${HIPMER_HUMAN_DATA}/s_1_2_sequence.permuted ]
    [ -f ${HIPMER_HUMAN_DATA}/fos1.fastq ] && [ -f ${HIPMER_HUMAN_DATA}/fos2.fastq ]
    [ -f ${HIPMER_HUMAN_DATA}/jmp1.fastq ] && [ -f ${HIPMER_HUMAN_DATA}/jmp2.fastq ]
fi
for i in ${HIPMER_HUMAN_DATA}/* ; do ln -f $i $RUNDIR || ln -fs $i $RUNDIR ; done

cp -p ${INST}/etc/meraculous/pipeline/meraculous-human*.config $RUNDIR 



