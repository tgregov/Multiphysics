#!/bin/sh
#Launch this in MinGW64 shell from MSYS2
#You should have installed the following package from MSYS2:
#mingw-w64-x86_64-cmake
#mingw-w64-x86_64-toolchain

export GMSHSDK=/c/tools/gmsh-4.6.0-Windows64-sdk #put gmsh sdk here

export PATH=${GMSHSDK}/bin:${GMSHSDK}/lin:${PATH}
export INCLUDE=${GMSHSDK}/include
export LIB=${GMSHSDK}/lib

cd ../../

rm -rf build
mkdir build
cd build

cmake -G "CodeBlocks - MinGW Makefiles" -DCMAKE_SH=SH-NOTFOUND  ..
cp -r ${GMSHSDK}/lib/gmsh-4.6.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgcc_s_seh-1.dll bin
cp -r /c/tools/msys64/mingw64/bin/libgomp-1.dll bin
cp -r /c/tools/msys64/mingw64/bin/libstdc++-6.dll bin
cp -r /c/tools/msys64/mingw64/bin/libwinpthread-1.dll bin

cd ../