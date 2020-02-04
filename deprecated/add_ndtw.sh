#!/bin/bash

source ${HOME}/.oh-my-zsh/custom/functions.zsh
load_spack

# Make sure SPACK_ROOT (spack root directory) is an environment variable
if [ -z ${SPACK_ROOT+x} ]; then
    echo "Please set SPACK_ROOT variable. Exiting configure script.";
    exit 1
fi

cwd=`pwd`
mpibin=`spack location -i openmpi`/bin/mpiexec

declare -a fdirs=("103x28" "205x55" "409x109" "817x217" "1633x433")

for fdir in "${fdirs[@]}"
do
cd "$fdir"
echo `pwd`
${mpibin} -np 1 ${HOME}/wind/NaluWindUtils/build/src/preprocessing/nalu_preprocess -i add_ndtw.yaml
cd "$cwd"
done
