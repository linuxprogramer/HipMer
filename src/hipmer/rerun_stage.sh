#!/bin/bash

if [ "$1" == "" ]; then
	echo "rerun_stage.sh <comma separated stages>"
	echo "  available stages are"
	for log in `grep "#*.log" timings.log|awk '{print $2}'`; do
		echo "    ${log%.log}"
	done
	echo "    last"
	exit 0
fi 

install_path=`grep "# Install path" timings.log|cut -f2 -d':'|xargs`

arrstage=(${1//,/ })

for s in ${arrstage[@]}; do 
    echo "Removing $s.log..."
    if [ "$s" != "last" ]; then
        rm -f "$s.log"
    fi
done


cfile=`ls meraculous-*.config`
export RUNDIR="."
export TIMINGS="timings-rerun.log"
$install_path/bin/run_hipmer.sh $install_path $cfile |grep -v "Skipping"
unset TIMINGS
unset RUNDIR

