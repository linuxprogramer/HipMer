#!/bin/bash -l

#module rm PrgEnv-intel
#module load PrgEnv-gnu
module load darshan
module load cmake
module load perl
#module load gnuplot
module load boost
module load bupc-narrow
# ignore failure here because it happens during builds 
module load lustre-cray_ari_s 2>/dev/null

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=edison
export BUILD=${BUILD:=$TMPDIR/build-${USER}-${HIPMER_BUILD_ENV}}
#export BUILD=${BUILD:=$SCRATCH/build-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
export CC=$(which cc)
export CXX=$(which CC)
export HIPMER_BUILD_OPTS="-DHIPMER_KMER_LENGTHS='19;21;31;51;61;71'"
export BUILD_TYPE="Release"
export CORES_PER_NODE=${CORES_PER_NODE:=24}

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"

