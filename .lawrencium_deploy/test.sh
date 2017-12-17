#!/bin/bash -l
[ -n "${HIPMER_BUILD_ENV}" ] || . $(dirname $0)/env.sh

set -ex

[ -n "${BUILD}" ]

cd ${BUILD}
make test




