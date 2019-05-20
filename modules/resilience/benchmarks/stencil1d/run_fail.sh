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
NTILES=128
#NTILES=32
NST=128
#NST=32
NS=1048576
#NS=65536

if [ "$#" -ge 1 ]; then
    NUM=$1
fi
if [ "$#" -ge 2 ]; then
    if [ "$2" -eq 1 ]; then
        NTILES=128
        NST=128
        NS=1048576
    elif [ "$2" -eq 0 ]; then
        NTILES=32
        NST=32
        NS=65536
    else
        echo "Unknown input size"
        exit
    fi
fi

echo "Times " $NUM
echo "Size " $NTILES $NS

echo "Run with checksum"

for i in $(seq $NUM);
do
    echo $i
    ./lw_1d_full_tile_replay -ntiles $NTILES -tilesize 16000 -nstepstask $NST -nsteps $NS -checksum
done

for failure_files in failure_files/replay/test_32_128-1-1_8192-its_0_f0.01_Work-stealing_balanced_lifo.out_failedTasks.txt failure_files/replay/test_32_128-1-1_8192-its_0_f0.1_Work-stealing_balanced_lifo.out_failedTasks.txt
do
    for i in $(seq $NUM);
    do
        echo $i
        ./lw_1d_full_tile_replay -ntiles $NTILES -tilesize 16000 -nstepstask $NST -nsteps $NS -checksum -inject $failure_files
    done
done


for i in $(seq $NUM);
do
    echo $i
    ./lw_1d_full_tile_replic -ntiles $NTILES -tilesize 16000 -nstepstask $NST -nsteps $NS -checksum
done

for failure_files in failure_files/replication/test_32_128-1-1_8192-its_0_f0.01_Work-stealing_balanced_lifo.out_failedTasks.txt failure_files/replication/test_32_128-1-1_8192-its_0_f0.1_Work-stealing_balanced_lifo.out_failedTasks.txt
do
    for i in $(seq $NUM);
    do
        echo $i
        ./lw_1d_full_tile_replic -ntiles $NTILES -tilesize 16000 -nstepstask $NST -nsteps $NS -checksum -inject $failure_files
    done
done

