#!/bin/sh
if [ ! -d "gmsh-4.1.4-Linux64-sdk" ]; then
  wget http://gmsh.info/bin/Linux/gmsh-4.1.4-Linux64-sdk.tgz
  tar -xf gmsh-4.1.4-Linux64-sdk.tgz 
  rm -rf gmsh-4.1.4-Linux64-sdk.tgz 
fi

cd gmsh-4.1.4-Linux64-sdk/
module load cmake/3.11.1
module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH}  

cd ../
rm -rf build/  
mkdir build

cd build/
cmake ../ -DCMAKE_BUILD_TYPE=Release  -G "Unix Makefiles" 
make
