#!/bin/sh
# Test on mac
clear
cd ../..
./build/bin/main ./Geometry/2D\ Rectangle/rectangle.msh ./Params/param.dat
#cd ./run/macOS/
gmsh results.msh
cd ./run/macOS