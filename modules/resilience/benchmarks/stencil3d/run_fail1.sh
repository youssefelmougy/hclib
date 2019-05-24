#!/bin/bash

NUM=1
NTILES=128
#NTILES=32
NST=128
#NST=32
NS=1048576
#NS=65536

echo "Times " $NUM

echo "Run no res no checksum"
for i in $(seq $NUM);
do
    ./heat_3d_simple_replay -ntiles 16 -tilesize 32 -nsteps 1024 -nochecksum
done 

echo "Run with checksum"
for failure_file in failure_files/replay/test_32_16-16-16_1024-its_0_f0.001_Work-stealing_balanced_lifo.out_failedTasks.txt failure_files/replay/test_32_16-16-16_1024-its_0_f0.01_Work-stealing_balanced_lifo.out_failedTasks.txt failure_files/replay/test_32_16-16-16_1024-its_0_f0.1_Work-stealing_balanced_lifo.out_failedTasks.txt  
do
    for i in $(seq $NUM);
    do
        ./heat_3d_simple_replay -ntiles 16 -tilesize 32 -nsteps 1024 -checksum -inject $failure_file
    done
done

