#!/bin/bash

NUM=1
SIZE=12000
TILE=400

if [ "$#" -ge 1 ]; then
    NUM=$1
fi
if [ "$#" -ge 2 ]; then
    if [ "$2" -eq 1 ]; then
        SIZE=24000
    elif [ "$2" -eq 0 ]; then
        SIZE=12000
    else
        echo "Unknown input size"
        exit
    fi
fi

echo "Times " $NUM
echo "Size " $SIZE
echo "Tile " $TILE

echo "Run with checksum"
for name in resilient_cholesky_replay resilient_cholesky_replic resilient_cholesky_abft
do
    for i in $(seq $NUM);
    do
        echo $i
        ./$name $SIZE $TILE
    done
    for rate in  0.01 0.1
    do
        for i in $(seq $NUM);
        do
            echo $i
            ./$name $SIZE $TILE $rate
        done
    done
done
