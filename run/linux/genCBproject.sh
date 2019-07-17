#!/bin/sh
#Do not forget to install libgfortran3  on your system!
HERE=$PWD

cd $HOME
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

cd $HERE
export GMSHSDK=${HOME}/gmsh-4.4.0-Linux64-sdk/
export EIGENSDK=${HOME}/eigen-eigen-323c052e1731/

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lib:${PATH}
export INCLUDE=${GMSHSDK}/include:${INCLUDE}
export INCLUDE=${EIGENSDK}:${INCLUDE}
export LIB=${GMSHSDK}/lib:${LIB}
export PYTHONPATH=${GMSHSDK}/lib:${PYTHONPATH}
export DYLD_LIBRARY_PATH=${GMSHSDK}/lib:${DYLD_LIBRARY_PATH}

cd ../../

rm -rf build
mkdir build
cd build

mkdir Release
cd Release
cmake ../../ -DCMAKE_BUILD_TYPE=Release  -G "CodeBlocks - Unix Makefiles"
cp -r ../../geometry/ $PWD/bin
cp -r ../../params/ $PWD/bin

cd ../

mkdir Debug
cd Debug
cmake ../../ -DCMAKE_BUILD_TYPE=Debug  -G "CodeBlocks - Unix Makefiles"
cp -r ../../geometry/ $PWD/bin
cp -r ../../params/ $PWD/bin

cd ../../
