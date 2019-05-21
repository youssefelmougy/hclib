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

echo "Run no res"
for i in $(seq $NUM);
do
    echo $i
    ./cg_baseline crankseg_1/crankseg_1.mtx 128 500 first_eig_value_vec.txt
done 

echo "Run with res"
for name in cg_replay cg_replication
do
    for i in $(seq $NUM);
    do
        echo $i
        ./$name crankseg_1/crankseg_1.mtx 128 500 failure_files/replay/test_32_spmv_crankseg_1_128-tiles_500-its_0_f0.0_Work-Stealing_lifo.out_failedTasks.txt first_eig_value_vec.txt 
    done
done

