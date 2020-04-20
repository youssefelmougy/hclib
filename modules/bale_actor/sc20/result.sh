#!/bin/bash

max_node=4
ranks=2
numtimes=5

while getopts i:n:r:m:h option
do
  case "${option}"
  in
    i) FILE=${OPTARG};;
    n) numtimes=${OPTARG};;
    r) ranks=${OPTARG};;
    m) max_node=${OPTARG};;
    h) echo "Options
             -i <input file>
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -r <1..inf> default is 2 i.e. two ranks per node
             -m <1..inf> default is 4 i.e. max number of nodes is 4
             -h show options"; exit 0;;
  esac
done

EXES=(histo ig permute randperm toposort transpose triangle)

count=$(cat $FILE |grep seconds|wc -l)
#echo $count

num_node_count=$(echo $max_node | awk '{print log($1)/log(2)}')
actual_count=$((numtimes*num_node_count*3*7))
#echo $actual_count

if [ $count -ne $actual_count ]; then
echo
echo "WARNING: Requires " $actual_count " readings but only " $count " found"
echo
fi

avg_times=$(cat $FILE |grep seconds|awk -v n_times="$numtimes" '{sum+=$1} NR%n_times==0 {print sum/n_times; sum=0}')
#echo $avg_times

i=0
bm_count=$((3*num_node_count))
bm_id=0
node_count=2

echo
echo "            HClib Actor/Selector"
echo
echo "        AGI        Conveyor     Selector"
for str in $avg_times
do
    #print benchmark name
    rem=$((i%bm_count))
    if [ "$rem" = 0 ]; then
    echo
        echo ${EXES[$bm_id]}
        bm_id=$((bm_id+1))
        node_count=2
    fi

    #print time
    printf "%12.4f  " $str

    #print node count after 3 reading
    i=$((i+1))
    rem=$((i%3))
    if [ "$rem" = 0 ]; then
        num_ranks=$((node_count*ranks))
        printf "%4d PEs\n" $num_ranks
        node_count=$((node_count*2))
    fi
done

echo

