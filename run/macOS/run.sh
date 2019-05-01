#!/bin/sh
# Test on mac
clear
cd ../..
./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_5.msh ./Params/param.dat ./results/results0_5.msh
#cd ./run/macOS/
#gmsh results/results0_1.msh
cd ./run/macOS