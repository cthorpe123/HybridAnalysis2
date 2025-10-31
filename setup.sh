#!/bin/bash

default_ubc='v10_04_07_15'

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh 

if [ -z "${1}" ]; then
    ubc=$default_ubc
else 
    ubc=$1
fi

echo "Setting up uboonecode ${ubc}"

setup uboonecode $ubc -q e26:prof
unsetup mrb
setup mrb -o

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$PWD/Funcs/
