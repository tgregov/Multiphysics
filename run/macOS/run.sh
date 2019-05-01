#!/bin/sh
# Test on mac
clear
cd ../..
./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_5.geo ./Geometry/2D\ Rectangle/newMesh.msh ./Params/param.dat ./results/newResults.msh
#cd ./run/macOS/
#gmsh results/newResults.msh
cd ./run/macOS