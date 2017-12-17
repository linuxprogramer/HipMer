#!/bin/bash -l

module rm PrgEnv-intel
module load PrgEnv-cray
module load darshan
module load gcc/4.8.1
module load cmake
module load cray-hdf5-parallel
module load perl
module load gnuplot
module load boost/1.55
module load lustre-cray_ari_s

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=edison-cray-upc
export BUILD=${BUILD:=$TMPDIR/build-${USER}-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
export CC=$(which gcc)
export CXX=$(which g++)
export CORES_PER_NODE=${CORES_PER_NODE:=24}
export HIPMER_BUILD_OPTS="-DCMAKE_UPC_COMPILER_INIT='$(which cc) -h upc'"

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"

