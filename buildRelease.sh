#!/bin/bash
set -e

if [ -d .git ]
then
  HIPMER_VERSION=$(git describe --tag)
  echo ${HIPMER_VERSION} >HIPMER_VERSION
fi
HIPMER_VERSION=$(cat HIPMER_VERSION)
MYCODE=$(pwd)

cd $TMPDIR
buildname=HipMer-${HIPMER_VERSION}
echo "Building ${buildname}"
rm -rf ${buildname}

set -x
git clone ${MYCODE} ${buildname}
cp -p ${MYCODE}/HIPMER_VERSION ${buildname}
rm -r ${buildname}/.git
tar=${MYCODE}/opt/${buildname}.tar.gz
tar -czf ${tar} ${buildname}
set +x

rm -rf ${buildname}
echo "Release is in ${tar}"

