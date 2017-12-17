#!/bin/bash -l

module rm PrgEnv-intel
module rm PrgEnv-gnu
module rm PrgEnv-cray
module load PrgEnv-cray
#module load git
module load darshan
module load cmake
module load perl
#module load gnuplot
module load lustre-cray_ari_s

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=cori-cray-upc
export BUILD=${BUILD:=/tmp/build-${USER}-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
export CC=$(which cc)
export CXX=$(which CC)
export CORES_PER_NODE=${CORES_PER_NODE:=32}
export HIPMER_BUILD_OPTS="-DCMAKE_UPC_COMPILER_INIT='$(which cc)\ -h\ upc' -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX"

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"

