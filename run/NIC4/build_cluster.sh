#!/bin/sh
if [ ! -d "gmsh-4.4.0-Linux64-sdk" ]; then
  wget http://gmsh.info/bin/Linux/gmsh-4.4.0-Linux64-sdk.tgz
  tar -xf gmsh-4.4.0-Linux64-sdk.tgz 
  rm -rf gmsh-4.4.0-Linux64-sdk.tgz 
fi

if [ ! -d "eigen-eigen-323c052e1731" ]; then
  wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
  tar -xf 3.3.7.tar.gz
  rm -rf 3.3.7.tar.gz
fi

if [ $1 == "vega" ]; then
  module load GCC/4.9.2 
elif [ $1 == "nic4" ]; then
  module load cmake/3.11.1
  module load gcc/4.9.2
fi
export CC=gcc
export CXX=g++

cd gmsh-4.4.0-Linux64-sdk/
export FC=gfortran
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH} 
cd ../

cd eigen-eigen-323c052e1731/
export INCLUDE=${PWD}:${INCLUDE}

cd ../
rm -rf build/  
mkdir build

cd build/
cmake ../ -DCMAKE_BUILD_TYPE=Release  -G "Unix Makefiles" 
make
