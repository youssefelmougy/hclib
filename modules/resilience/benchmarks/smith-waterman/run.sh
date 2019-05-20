#!/bin/bash
#SBATCH -q debug
##SBATCH -q regular
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:30:00
##SBATCH -t 00:98:00

#export HCLIB_LOCALITY_FILE=$HOME/resilience/amt-resilience/src/apps/examples/sandia.2socket.json
#export HCLIB_WORKERS=32

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

for name in smith_waterman.baseline.out  smith_waterman.replay.out  smith_waterman.replication.out
do
    for i in $(seq $NUM);
    do
        echo $i
        sh run_base.sh $SIZE $name
    done
done

