#!/bin/bash -l

module purge

module load uge
module load PrgEnv-gnu/4.8
module load cmake
module load OFED
module load openmpi/1.6.5
module load perl
module load gnuplot
module load boost/1.57.0
module load zlib

module use ~regan/modulefiles-genepool
module load bupc/regan-2.22.0-gnu48-sysv

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi


export HIPMER_BUILD_ENV=genepool-debug
export BUILD_TYPE=Debug
export BUILD=${BUILD:=/tmp/build-${USER}-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
#export HIPMER_BUILD_OPTS="-DCMAKE_UPC_USE_PTHREADS=8 -DHIPMER_VERBOSE=2 "
export HIPMER_BUILD_OPTS="-DHIPMER_VERBOSE=1"

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"
