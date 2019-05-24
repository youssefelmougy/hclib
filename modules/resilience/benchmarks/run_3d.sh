#!/bin/bash

BIG=1;
NUM=5;

while getopts b:n:h option
do
  case "${option}"
  in
    b) BIG=${OPTARG};;
    n) NUM=${OPTARG};;
    h) echo "Options
             -b <0|1> default is 1 i.e use bigger inputs (same as paper) which might take more time
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -h show options"; exit 0;;
  esac
done

echo use big inputs $BIG
echo number of times $NUM

pushd ../../../
./clean.sh
./install.sh
source hclib-install/bin/hclib_setup_env.sh
popd

cd stencil3d
make clean
make
echo "## start stencil3d ##"
/bin/bash -x run.sh $NUM $BIG
echo "## end stencil3d ##"
echo
echo "## start stencil3d mix ##"
/bin/bash -x run_mix.sh $NUM $BIG
echo "## end stencil3d mix ##"
echo
echo "## start stencil3d fail ##"
/bin/bash -x run_fail.sh $NUM $BIG
echo "## end stencil3d fail ##"
cd ..

pushd ../../../
./clean.sh
popd

