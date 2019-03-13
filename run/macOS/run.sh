#!/bin/sh
# Test on mac
clear
cd ../..
./build/bin/main ./Geometry/2D\ Rectangle/rectangle.msh
#cd ./run/macOS/
gmsh results.msh