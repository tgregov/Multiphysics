#!/bin/sh
wget http://gmsh.info/bin/Linux/gmsh-4.4.0-Linux64-sdk.tgz
tar -xf gmsh-4.4.0-Linux64-sdk.tgz 
cd gmsh-4.4.0-Linux64-sdk/
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH}
cd ../
wget http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz
tar -xf 3.3.7.tar.gz
cd ./eigen-eigen-323c052e1731
export INCLUDE=${PWD}:${INCLUDE}
cd ../
git clone https://github.com/tgregov/Multiphysics.git
cd Multiphysics/  
mkdir build
cd build 
cmake ../ -DCMAKE_BUILD_TYPE=Release  -G "Unix Makefiles" 
make
