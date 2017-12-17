#!/bin/bash -l

module purge

module load uge
module load PrgEnv-gnu/4.8
module load cmake
module load hdf5
module load perl
module load gnuplot
module load boost/1.57.0
module load zlib

module list

if [ -n "${HIPMER_BUILD_ENV}" ]
then
  unset BUILD
  unset PREFIX
fi


export HIPMER_BUILD_ENV=genepool-no-mpi
export BUILD=${BUILD:=${SCRATCH}/build-${USER}-${HIPMER_BUILD_ENV}}
export PREFIX=${PREFIX:=$SCRATCH/install-${HIPMER_BUILD_ENV}}
export HIPMER_INSTALL=${HIPMER_INSTALL:=${PREFIX}}

export HIPMER_POST_INSTALL="cp ${BASH_SOURCE[0]} ${PREFIX}"
