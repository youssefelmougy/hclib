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

echo "Run with checksum"

for i in $(seq $NUM);
do
    echo $i
    ./heat_3d_simple_replay -ntiles $NTILES -tilesize $NST -nsteps $NS -checksum
done

for failure_file in failure_files/replay/test_32_16-16-16_1024-its_0_f0.01_Work-stealing_balanced_lifo.out_failedTasks.txt failure_files/replay/test_32_16-16-16_1024-its_0_f0.1_Work-stealing_balanced_lifo.out_failedTasks.txt  
do
    for i in $(seq $NUM);
    do
        echo $i
        ./heat_3d_simple_replay -ntiles $NTILES -tilesize $NST -nsteps $NS -checksum -inject $failure_file
    done
done


for i in $(seq $NUM);
do
    echo $i
    ./heat_3d_simple_replic -ntiles $NTILES -tilesize $NST -nsteps $NS -checksum
done

for failure_file in failure_files/replication/test_32_16-16-16_1024-its_0_f0.01_Work-stealing_balanced_lifo.out_failedTasks.txt failure_files/replication/test_32_16-16-16_1024-its_0_f0.1_Work-stealing_balanced_lifo.out_failedTasks.txt
do
    for i in $(seq $NUM);
    do
        echo $i
        ./heat_3d_simple_replic -ntiles $NTILES -tilesize $NST -nsteps $NS -checksum -inject $failure_file
    done
done

