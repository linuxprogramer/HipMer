#grep log timings.log|sort |awk '{print $2, $3}'> /dev/shm/lala
#grep Overall *.log*|sort|awk '{print $6}'> /dev/shm/gaga

fname="/tmp/$USER-stage-timings"

printf "%-35s %10s %10s %10s\n" "STAGE" "EXTERNAL" "INTERNAL" "DIFFERENCE" 

logs=`grep log timings.log|awk '{print $2}'` 
rm -f $fname
onoSetId=1
round=1
nstages=0
for log in $logs; do 
    if [[ $log == per_thread/oNo* ]]
    then
        echo -n "-->" >> $fname
    else
        nstages=$((nstages+1))
    fi
    #grep "# $log " timings.log|awk '{printf "%s %s", $2, $3}' >> $fname
    echo -n `basename $log |cut -d'.' -f1` >> $fname
    grep "# $log " timings.log|awk '{printf " %s", $3}' >> $fname
    grep Overall $log|awk '{printf " %s\n", $6}' >> $fname
    if [ "$log" == "per_thread/merger-links-${onoSetId}.log" ]
    then
        # get the max internal value
        internal=`grep Overall oNo-${onoSetId}-* | awk '{print $6}' |sort -nr |head -n 1`
        external=`grep "# parallel_oNo-${onoSetId}" timings.log|awk '{print $3}'`
        difference=`echo "$external-$internal"|bc -l`
        echo parallel_oNo-$onoSetId $external $internal $difference >> $fname
        onoSetId=$((onoSetId+1))
    fi
done

a1=`egrep -v "^-->" $fname |awk '{sum+=$2} END {if (NR>0) print sum}'`
a2=`egrep -v "^-->" $fname |awk '{sum+=$3} END {if (NR>0) print sum}'`
echo "TOTAL $a1 $a2" >> $fname

cat $fname |awk '{printf "%-35s %10s %10s %10s\n", $1, $2, $3, $2-$3}' 

avdiff=`echo "($a1-$a2)/$nstages" | bc -l`
printf "Number of stages: %d ; average difference per stage: %.0f\n" $nstages $avdiff

tail -n1 timings.log

rm -f $fname

