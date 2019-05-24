#!/bin/bash

NUM=1
NTILES=16
#NTILES=16
NST=32
#NST=16
NS=1024
#NS=128

if [ "$#" -ge 1 ]; then
    NUM=$1
fi
if [ "$#" -ge 2 ]; then
    if [ "$2" -eq 1 ]; then
        NTILES=16
        NST=32
        NS=1024
    elif [ "$2" -eq 0 ]; then
        NTILES=16
        NST=16
        NS=128
    else
        echo "Unknown input size"
        exit
    fi
fi

echo "Times " $NUM
echo "Size " $NTILES $NS

echo "Run mix with checksum"
for ratio in 0 20 40 60 80 100
do
for i in $(seq $NUM);
do
    echo $i
    ./heat_3d_simple_mix -ntiles $NTILES -tilesize $NST -nsteps $NS -checksum -replic $ratio
done 
done

