#!/bin/bash -l


module rm PrgEnv-pgi
module load PrgEnv-intel
module load cmake
module load perl
module load gnuplot
module load bupc-narrow

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=hopper-debug
export BUILD=${BUILD:=$SCRATCH/build-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
export CC=$(which cc)
export CXX=$(which CC)
export CORES_PER_NODE=${CORES_PER_NODE:=24}
export HIPMER_BUILD_OPTS="-DHIPMER_ONLY=1"
export BUILD_TYPE="Debug"

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"

