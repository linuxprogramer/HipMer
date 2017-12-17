#!/bin/bash -l
set -e

[ -n "${HIPMER_BUILD_ENV}" ] || . $(dirname $0)/env.sh

[ -n "${BUILD}" ]
[ -n "${PREFIX}" ]

SRC=$(pwd)
if [ -e "${BUILD}" ]
then 
    cmakecache=${BUILD}/CMakeCache.txt
    if [ -f ${cmakecache} ]
    then
      testsame=$( ( grep HipMer_SOURCE_DIR /tmp/build-regan-edison/CMakeCache.txt ; echo "HipMer_SOURCE_DIR:STATIC=${SRC}" ) | uniq | wc -l)
      if [ "${testsame}" != "1" ]
      then
        echo "Source dirs do not match.  performing a DIST_CLEAN build"
        DIST_CLEAN=1
      fi
    fi

    if [ -n "${DIST_CLEAN}" ]
    then
      chmod -R u+w ${BUILD}
      rm -r ${BUILD}
      mkdir ${BUILD}
      rm -rf ${PREFIX}
    elif [ -n "${CLEAN}" ]
    then
      (cd ${BUILD} ; make clean )
    fi

else
    mkdir ${BUILD}
fi

BUILD_TYPE=${BUILD_TYPE:=Release}
cd ${BUILD}
if [ -x /usr/bin/lfs ]
then
  /usr/bin/lfs setstripe -c 1 . || /bin/true
fi
export HDF5_ROOT=${HDF5_DIR}
cmakelog=${BUILD}/cmake.log
export TMPDIR=/tmp
if [ ! -f ${cmakelog} ]
then
  time cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${HIPMER_BUILD_OPTS} ${SRC} 2>&1 | tee -a ${cmakelog}.tmp \
     && mv ${cmakelog}.tmp ${cmakelog}
fi

make_threads=${BUILD_THREADS:=$(lscpu |grep "^CPU(s):"|awk '{print $2}')}
echo Using $make_threads threads for the build
time make -j ${make_threads} || TMPDIR=/tmp make VERBOSE=1 2>&1 | tee make.err




