#!/bin/bash
#SBATCH -q debug
##SBATCH -q regular
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:30:00
##SBATCH -t 00:58:00

#export HCLIB_LOCALITY_FILE=$HOME/resilience/amt-resilience/src/apps/examples/sandia.2socket.json
#export HCLIB_WORKERS=32

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

echo "Run no res no checksum"
for i in $(seq $NUM);
do
    echo $i
    ./heat_3d_simple_nonres -ntiles $NTILES -tilesize $NST -nsteps $NS -nochecksum
done 

echo "Run with checksum"
#for name in heat_3d_simple_nonres heat_3d_simple_replay heat_3d_simple_replic
for name in heat_3d_simple_replay heat_3d_simple_replic
do
    for i in $(seq $NUM);
    do
        echo $i
        ./$name -ntiles $NTILES -tilesize $NST -nsteps $NS -checksum
    done
done

