#!/bin/bash

NUM=5;

while getopts n:h option
do
  case "${option}"
  in
    n) NUM=${OPTARG};;
    h) echo "Options
             -n <1..inf> default is 5 i.e. run each experiment 5 times
             -h show options"; exit 0;;
  esac
done

echo number of times $NUM

echo "### Start Build ###"
pushd ../../../
./clean.sh
CFLAGS=-DMAX_NUM_WAITS=256 ./install.sh
source hclib-install/bin/hclib_setup_env.sh
popd
echo "###End Build ###"
echo

echo "### Start Experiment ###"
cd cg
make clean
make
rm -rf crankseg_1.tar.gz crankseg_1
wget https://sparse.tamu.edu/MM/GHS_psdef/crankseg_1.tar.gz
tar -xvzf crankseg_1.tar.gz
echo
echo "### Start no failure run ###"
/bin/bash -x run.sh $NUM
echo "### End no failure run ###"
echo
echo "### Start failure run ###"
/bin/bash -x run_fail.sh $NUM
echo "### End failure run ###"
cd ..
echo "### End Experiment ###"
echo

echo "### Start Cleanup ###"
pushd ../../../
./clean.sh
popd
echo "### End Cleanup ###"
echo

