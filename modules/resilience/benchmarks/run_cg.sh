#!/bin/bash
#SBATCH -q debug
##SBATCH -q regular
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:30:00
##SBATCH -t 00:58:00

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

pushd ../../../
./clean.sh
CFLAGS=-DMAX_NUM_WAITS=256 ./install.sh
source hclib-install/bin/hclib_setup_env.sh
popd

cd cg
make clean
make
rm -rf crankseg_1.tar.gz crankseg_1
wget https://sparse.tamu.edu/MM/GHS_psdef/crankseg_1.tar.gz
tar -xvzf crankseg_1.tar.gz
sh run.sh $NUM
sh run_fail.sh $NUM
cd ..

