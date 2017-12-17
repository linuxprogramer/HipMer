#!/bin/bash -l


module rm PrgEnv-intel
module load PrgEnv-gnu
module load darshan
module load cmake
module load perl
module load gnuplot
module load boost
module load bupc-narrow
module load lustre-cray_ari_s

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=edison-gnu
export BUILD=${BUILD:=$TMPDIR/build-${USER}-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
export CC=$(which cc)
export CXX=$(which CC)
export HIPMER_BUILD_OPTS="-DMPI_C_COMPILER=${CC} -DMPI_CXX_COMPILER=${CXX}"
export CORES_PER_NODE=${CORES_PER_NODE:=24}

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"

