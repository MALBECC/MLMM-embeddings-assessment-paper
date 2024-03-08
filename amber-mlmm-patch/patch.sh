#!/bin/bash

sander_path=$1
cp CMakeLists.txt ${sander_path}/.
cp qm2_extern_mlmm_module.F90 ${sander_path}/.
cp qm2_extern_module.F90 ${sander_path}/.

echo "The following files have been replaced:"
echo "${sander_path}/qm2_extern_module.F90"
echo "${sander_path}/CMakeLists.txt"
echo "The following file have been added to the sander path:"
echo "${sander_path}/qm2_extern_mlmm_module.F90

