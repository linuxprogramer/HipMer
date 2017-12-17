#!/bin/bash -l
basepath=$(dirname $0)

pids=
recurse_pkill()
{
  local killpids=($*)
  if [ ${#killpids[*]} -eq 0 ]
  then
    return
  fi
  echo "Finding Children of ${killpids[*]}"
  for pid in ${killpids[*]}
  do
    local children=($(pgrep -P ${pid}))
    if [ ${#children[*]} -gt 0 ]
    then
      recurse_pkill ${children[*]}
    fi
  done
  echo Killing ${killpids[*]}
  kill ${killpids[*]}
}
cleanup()
{
  set +e
  if [ -n "${pids}" ]
  then
    recurse_pkill ${pids}
  fi
}

BUILD_THREADS=${BUILD_THREADS:=$(lscpu |grep "^CPU(s):"|awk '{print $2}')}
builds=($(echo ${basepath}/env*.sh))
numbuilds=${#builds[@]}
if [ ${numbuilds} -gt 1 ]
then
  BUILD_THREADS=$((BUILD_THREADS / numbuilds + 1))
fi
export BUILD_THREADS

trap cleanup 0 1 2 3 15 
for env in ${builds[@]}
do
  ( (
    . $env
    ${basepath}/build.sh \
    && ${basepath}/install.sh \
    && ${basepath}/test.sh \
  ) > ${env##*/}-$(uname -n).log 2>&1 && echo "$env Built" ) || echo "$env Failed to build" &
  pid=$!
  echo "Started build for $env $pid"
  pids="${pids} ${pid}"
done

echo "Waiting for builds $(date): ${pids}"
wait ${pids}
trap '' 0
echo "Builds finished ($?) $(date)"


