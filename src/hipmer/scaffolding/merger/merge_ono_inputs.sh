#!/bin/bash

# We need to be able to run k separate instances of oNo simulataneously, so each oNo instance must 
# have N/k threads, where N is the number of cores in the run 
# Typically, k will be 9 (2 to 10)

round()
{
    echo $((($1+$2-1)/$2))
}

start_t=`date +%s.%N`
n_cores=$1
k=$2
threads_per_ono=$((n_cores/k))
fname=$3
if [ "$fname == meta" ]; then
    links_meta=$4
    if [ -n "$links_meta" ]; then 
        fname=`awk '{print $3}' $links_meta`
        insert=`awk '{print $1}' $links_meta`
        std=`awk '{print $2}' $links_meta`
        echo -e "$insert\t$std\tmerged_${fname}" > merged_$links_meta
        echo "New linksMeta: merged_$links_meta"
    fi
fi

echo "Merging files $fname from $n_cores threads into $threads_per_ono threads for $k oNo instances"

thread=0
cores=""
echo -n > merged_${fname}_$thread
for i in `seq 0 $((n_cores-1))`; do 
    cores="$cores $i"
    cat ${fname}_$i >> merged_${fname}_$thread
    #if [ $i != 0 ] && [ $((($i+1)%$k)) == 0 ]; then
    if [ $((($i+1)%$k)) == 0 ]; then
        if [ $thread -lt $((threads_per_ono-1)) ]; then 
            echo "Cat $cores written to merged_${fname}_$thread"
            thread=$((thread+1))
            cores=""
        fi
    fi
done

echo "Cat $cores written to merged_${fname}_$thread"
        
stop_t=`date +%s.%2N`
diff_t=`echo $stop_t-$start_t | bc -l`
cmd=`basename $0`
printf "Overall time for %s is %.2f s\n" $cmd $diff_t
