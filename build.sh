#!/bin/sh
wget http://gmsh.info/bin/Linux/gmsh-4.1.4-Linux64-sdk.tgz
tar -xf gmsh-4.1.4-Linux64-sdk.tgz 
cd gmsh-4.1.4-Linux64-sdk/include
rm -rf gmsh.h
mv gmsh.h_cwrap gmsh.h
cd ../ 
export PATH=${PWD}/bin:${PWD}/lib:${PATH}
export INCLUDE=${PWD}/include:${INCLUDE}
export LIB=${PWD}/lib:${LIB}
export PYTHONPATH=${PWD}/lib:${PYTHONPATH}
cd ../
git clone https://github.com/ImperatorS79/Multiphysics.git
cd Multiphysics/  
mkdir build
cd build 
cmake ../ -DCMAKE_BUILD_TYPE=Release  -G "Unix Makefiles" 
make
