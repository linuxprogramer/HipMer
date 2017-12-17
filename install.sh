#!/bin/bash


if [ $# -lt 1 ]
then
  echo "Usage: $0 <destination> " >&2
  exit 0
fi
set -e

installdir=$1
[ -e "${installdir}" ] || mkdir ${installdir}

CXX=g++
CC=gcc
arch=$(uname -s)
if [ "$arch" == "Darwin" ]
then

  CXX=clang++
  CC=clang
  if which g++-4.8 && which gcc-4.8 && [ -d /usr/local/include/boost ]
  then
    echo "Found g++-4.8 gcc-4.8 and a boost library"
    CXX=g++-4.8
    CC=gcc-4.8
  else
    echo "For Mac OS X: Expecting HomeBrew installed packages: argp-standalone boost gcc48 coreutils open-mpi"
    echo "brew install gcc48 argp-standalone boost hdf5 coreutils"
    echo "brew install --cc=gcc-4.8 --build-from-source open-mpi"
    echo "If build fails, please install these and try again"
  fi
  installdir=$(greadlink -f ${installdir})

else
  installdir=$(readlink -f ${installdir})
fi

if [ -d ./build ]
then
    echo "Error:"
    echo "Older build found in the current directory. Please delete/rename it before running the installer"
    exit 0
fi

mkdir build \
 && cd build \
 && CC=$CC CXX=$CXX cmake -DCMAKE_INSTALL_PREFIX=${installdir} -DCMAKE_BUILD_TYPE=Release .. \
 && (make -j4 || make)  \
 && make install \
 || (echo "ERROR: Could not build and install" ; exit 1 )

cd ..

