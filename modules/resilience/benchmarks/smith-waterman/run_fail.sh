#!/bin/bash

NUM=1
SIZE=huge

if [ "$#" -ge 1 ]; then
    NUM=$1
fi
if [ "$#" -ge 2 ]; then
    if [ "$2" -eq 1 ]; then
        SIZE=huge
    elif [ "$2" -eq 0 ]; then
        SIZE=larger
    else
        echo "Unknown input size"
        exit
    fi
fi

echo "Times " $NUM
echo "Size " $SIZE

for name in smith_waterman.replay.out  smith_waterman.replication.out
do
for rate in 0 0.01 0.1
do
    for i in $(seq $NUM);
    do
        echo $i
        /bin/bash -x run_fail_base.sh $SIZE $name $rate
    done
done
done

