#!/bin/sh
# Test on mac
clear
cd ../..
./build/bin/main ./Geometry/2D\ Rectangle/rectangleBlala.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/param.dat ./results/newResultsTest.msh
gmsh results/newResultsTest.msh
cd ./run/macOS