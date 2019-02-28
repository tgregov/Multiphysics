#!/bin/sh
# Compile on Mac
GMSHSDK=/Users/.../sdk 		# ADD YOUR PATH TO THE GMSH SDK
EIGENSDK=/Users/.../eigen 	# ADD YOUR PATH TO THE EIGEN LIBRARY
export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:${PATH}
export INCLUDE=${GMSHSDK}/include:${INCLUDE}
export INCLUDE=${EIGENSDK}:${INCLUDE}
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