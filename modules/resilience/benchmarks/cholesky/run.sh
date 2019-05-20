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

