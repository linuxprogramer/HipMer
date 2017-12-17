#!/bin/bash -l

module rm PrgEnv-intel
module rm PrgEnv-gnu
module rm PrgEnv-cray
module load PrgEnv-intel
module rm gcc
module load gcc/4.9.3
#module load git
module load darshan
module load cmake
module load perl
#module load gnuplot
module load bupc-narrow
module load lustre-cray_ari_s

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=cori
export BUILD=${BUILD:=/tmp/build-${USER}-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=$PREFIX}
export CC=$(which cc)
export CXX=$(which CC)
export BUILD_TYPE="Release"
export CORES_PER_NODE=${CORES_PER_NODE:=32}
export HIPMER_BUILD_OPTS="-DHIPMER_KMER_LENGTHS=19;21;31;51;61;71"

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"

