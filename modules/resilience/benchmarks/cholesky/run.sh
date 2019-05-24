#!/bin/bash

NUM=1

if [ "$#" -ge 1 ]; then
    NUM=$1
fi

echo "Times " $NUM

echo "Run with checksum"
for name in resilient_cholesky_nonres resilient_cholesky_replay resilient_cholesky_replic resilient_cholesky_abft
do
    for i in $(seq $NUM);
    do
        echo $i
        ./$name 12000 400
    done
done

