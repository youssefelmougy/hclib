#!/bin/bash

set -x

FIG=1;
BIG=1;
NUM=5;

while getopts f:b:n:h option
do
  case "${option}"
  in
    f) FIG=${OPTARG};;
    b) BIG=${OPTARG};;
    n) NUM=${OPTARG};;
    h) echo "Options
             -f <1|2|3|4> default is 1 i.e. generate figure 1
             -b <0|1> default is 1 i.e use bigger inputs (same as paper) which might take more time
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -h show options"; exit 0;;
  esac
done

echo figure $FIG
echo use big inputs $BIG
echo number of times $NUM

if [ $FIG -eq 1 ]; then
    echo "## start stencil1d ##"
    cd stencil1d
    /bin/bash run.sh $NUM $BIG
    cd ..
    echo "## end stencil1d ##"
    
    echo "## start smith-waterman ##"
    cd smith-waterman
    /bin/bash run.sh $NUM $BIG
    cd ..
    echo "## end smith-waterman ##"

    echo "## start cholesky ##"
    cd cholesky
    /bin/bash run.sh $NUM
    cd ..
    echo "## end cholesky ##"

fi

if [ $FIG -eq 2 ]; then
    echo "## start stencil1d ##"
    cd stencil1d
    /bin/bash run_mix.sh $NUM $BIG
    cd ..
    echo "## end stencil1d ##"

fi

if [ $FIG -eq 3 ]; then
    echo "## start stencil1d ##"
    cd stencil1d
    /bin/bash run_fail.sh $NUM $BIG
    cd ..
    echo "## end stencil1d ##"

    echo "## start smith-waterman ##"
    cd smith-waterman
    /bin/bash run_fail.sh $NUM $BIG
    cd ..
    echo "## end smith-waterman ##"

    echo "## start cholesky ##"
    cd cholesky
    /bin/bash run_fail.sh $NUM
    cd ..
    echo "## end cholesky ##"

fi
echo

