#!/bin/bash

NUM=1

if [ "$#" -ge 1 ]; then
    NUM=$1
fi

echo "Times " $NUM

echo "Run with res"
for name in replay replication
do
  for failure_file in failure_files/$name/test_32_spmv_crankseg_1_128-tiles_500-its_0_f0.0_Work-Stealing_lifo.out_failedTasks.txt failure_files/$name/test_32_spmv_crankseg_1_128-tiles_500-its_0_f0.01_Work-Stealing_lifo.out_failedTasks.txt failure_files/$name/test_32_spmv_crankseg_1_128-tiles_500-its_0_f0.1_Work-Stealing_lifo.out_failedTasks.txt
  do
    for i in $(seq $NUM);
    do
        echo $i
        ./cg_$name crankseg_1/crankseg_1.mtx 128 500 $failure_file first_eig_value_vec.txt 
    done
  done
done

