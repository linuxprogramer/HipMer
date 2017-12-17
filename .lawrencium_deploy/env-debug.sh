#!/bin/bash -l

module rm intel
module rm openmpi
module load intel/2016.1.150
module load openmpi
module load berkeley_upc/2.22.0-intel
module load perl
module load cmake
module load gnuplot

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi

export HIPMER_BUILD_ENV=lawrencium-debug
export BUILD=${BUILD:=$SCRATCH/build-${USER}-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}
export HIPMER_BUILD_OPTS="-DCMAKE_UPC_COMPILER_ENV_VAR=upcc -DMPI_C_COMPILER=$(which mpicc) -DMPI_CXX_COMPILER=$(which mpic++)"
export CORES_PER_NODE=${CORES_PER_NODE:=16}
export BUILD_TYPE="Debug"

export CC=$(which icc)
export CXX=$(which icpc)
export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"
