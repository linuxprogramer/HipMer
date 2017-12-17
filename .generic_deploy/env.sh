#!/bin/bash -l

[ -z "$SCRATCH" ] && echo "Define \$SCRATCH for the root of the build and install dirs" && exit 1

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=generic
#export BUILD=${BUILD:=$SCRATCH/build-${HIPMER_BUILD_ENV}}
#export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export BUILD=${SCRATCH}/build-${HIPMER_BUILD_ENV}
export PREFIX=${SCRATCH}/install-${HIPMER_BUILD_ENV}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
export CC=$(which cc)
export CXX=$(which CC)
export HIPMER_BUILD_OPTS=""
export BUILD_TYPE="Release"
#export CORES_PER_NODE=${CORES_PER_NODE:=4}
#export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"

