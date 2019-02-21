#!/bin/sh
# Compile on Mac
GMSHSDK=/Users/thomasgregov9/Documents/Unif/1Master\ Inge\ Civil/Cours/2Quadri/Multiphysics\ integrated\ computational\ project/Gmsh/sdk
export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:${PATH}
export INCLUDE=${GMSHSDK}/include:${INCLUDE}
export LIB=${GMSHSDK}/lib:${LIB}
export PYTHONPATH=${GMSHSDK}/lib:${PYTHONPATH}
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:${DYLD_LIBRARY_PATH}

cd ..
cd ..
rm -rf ./build
mkdir build
cd build
cmake ..
make
cd ..
cd ./run/macOS/