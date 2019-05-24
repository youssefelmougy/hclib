#!/bin/bash

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

MPIRUN=mpirun
#MPIRUN=srun

for nnodes in 2 4
do

    echo "Run no res no checksum"
    for i in $(seq $NUM);
    do
        echo $i
        $MPIRUN -n $nnodes ./lw_1d_full_tile_mpi_nonres -ntiles $(( $NTILES * $nnodes )) -tilesize 16000 -nstepstask $NST -nsteps $NS -nochecksum
    done 
    
    echo "Run with checksum"
    #for name in lw_1d_full_tile_mpi_nonres lw_1d_full_tile_mpi_replay lw_1d_full_tile_mpi_replic
    for name in lw_1d_full_tile_mpi_replay lw_1d_full_tile_mpi_replic
    do
        for i in $(seq $NUM);
        do
            echo $i
            $MPIRUN -n $nnodes ./$name -ntiles $(( $NTILES * $nnodes )) -tilesize 16000 -nstepstask $NST -nsteps $NS -checksum
        done
    done

done

