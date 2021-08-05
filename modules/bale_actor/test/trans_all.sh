#!/bin/bash

for app in histo ig permute randperm transpose triangle toposort;
do
    ./trans.sh ${app}
    make ${app}_selector_lambda.trans
done
