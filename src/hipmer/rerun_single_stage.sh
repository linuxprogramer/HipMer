#!/bin/bash

if [ "$1" == "meraculous" ]; then
	stage=contigs.log
elif [ "$1" == "last" ]; then
	stage=last
else 
	stage=per_thread/$1.log
fi

if [ "$stage" == "last" ]; then
	x=`tail -n1 timings.log`
    complete=`tail -n1 timings.log|cut -f2-2 -d' '`
    if [ "$complete" == "Total" ]; then 
        x=`tail -n4 timings.log`
    fi
	y=${x##'# Starting '}
	cmd=${y%at*}
else
	cmd=`grep -B1 " $stage" timings.log|cut -f1 -d'#'`
fi

if [ "$1" == "" ] || [ "$cmd" == "" ]; then
	echo "rerun_stage.sh <stage>"
	echo "  available stages are"
	for log in `grep "#*.log" timings.log|awk '{print $2}'`; do
        log=`basename $log`
		if [ "$log" == "contigs.log" ]; then
			echo "    meraculous"
		else
			echo "    ${log%.log}"
		fi
	done
	echo "    last"
	exit 0
fi 

#echo $cmd

install_path=`grep "# Install path" timings.log|cut -f2 -d':'|xargs`
#echo $install_path

prev_path=$PATH
export PATH=$PATH:$install_path/bin
read -e -p "" -i "$cmd" newcmd
$newcmd
export PATH=$prev_path

#$cmd
