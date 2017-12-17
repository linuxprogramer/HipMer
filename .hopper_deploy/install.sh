#!/bin/bash -l
set -x

[ -n "${HIPMER_BUILD_ENV}" ] || . $(dirname $0)/env.sh

set -ex

[ -n "${BUILD}" ]
[ -n "${PREFIX}" ]

[ -z "${CLEAN}" ] || rm -fr ${PREFIX}
mkdir -p ${PREFIX}
if [ -x /usr/bin/lfs ]
then
  /usr/bin/lfs setstripe -c 1 ${PREFIX} || true
fi
cd ${BUILD}
TMPDIR=/tmp make -j 16 install
cd -
deploy=${BASH_SOURCE[0]}
deploy=${deploy%/*}
[ ! -d ${deploy}/module_dependencies ] || cp -p ${deploy}/module_dependencies $PREFIX/.deps

${HIPMER_POST_INSTALL} || /bin/true
