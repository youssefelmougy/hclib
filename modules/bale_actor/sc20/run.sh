#!/bin/bash

max_node=4
ranks=2
size=s
numtimes=5

while getopts s:n:r:m:h option
do
  case "${option}"
  in
    s) size=${OPTARG};;
    n) numtimes=${OPTARG};;
    r) ranks=${OPTARG};;
    m) max_node=${OPTARG};;
    h) echo "Options
             -s <s|m|l>  default is s i.e. use small input
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -r <1..inf> default is 2 i.e. two ranks per node
             -m <1..inf> default is 4 i.e. max number of nodes is 4
             -h show options"; exit 0;;
  esac
done

EXES=(histo ig permute randperm toposort transpose triangle)

if [ "$SIZE" = s ]; then
    echo "Using small input"
    NS=(10000 10000 1000 10000 1000 1000 100)
elif [ "$SIZE" = m ]; then
    echo "Using medium input"
    NS=(1000000 1000000 10000 100000 10000 10000 1000)
elif [ "$SIZE" = l ]; then
    echo "Using large input"
    NS=(10000000 10000000 100000 1000000 100000 100000 10000)
else
    echo "Using small input (default)"
    NS=(10000 10000 1000 10000 1000 1000 100)
fi

for i in $(seq 0 6)
do

  EXE=${EXES[$i]}
  N=${NS[$i]}
  echo $EXE $N $i

  nodes=2
  while [ $nodes -le $max_node ]; do

    count=$((nodes*ranks))

    echo "At $i , nodes $nodes"

    echo "${EXE}_agi"

    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_agi $c"
      $OSHRUN -n $count ./${EXE}_agi -n $N
      sleep 10
    done

    echo "${EXE}_conveyor"
    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_conveyor $c"
      $OSHRUN -n $count ./${EXE}_conveyor -n $N
      sleep 10
    done

    echo "${EXE}_selector"

    for c in $(seq 1 $numtimes)
    do
      echo "in_${EXE}_selector $c"
      $OSHRUN -n $count ./${EXE}_selector -n $N
      sleep 10
    done

    nodes=$((nodes*2))
  done
done

