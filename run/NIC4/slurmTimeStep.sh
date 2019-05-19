#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=06:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=terminal.txt
##SBATCH --mail-user=joachim.marichal@student.uliege.be
#SBATCH --mail-type=ALL

module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=16

cd $HOME/Multiphysics

./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_025.geo ./Geometry/2D\ Rectangle/newMesh4.msh  ./Params/paramT4.dat ./results/newResults4.msh ./results/errorLastTest.txt
./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_025.geo ./Geometry/2D\ Rectangle/newMesh5.msh  ./Params/paramT5.dat ./results/newResults5.msh ./results/errorLastTest.txt
./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_025.geo ./Geometry/2D\ Rectangle/newMesh6.msh  ./Params/paramT6.dat ./results/newResults6.msh ./results/errorLastTest.txt
./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_025.geo ./Geometry/2D\ Rectangle/newMesh7.msh  ./Params/paramT7.dat ./results/newResults7.msh ./results/errorLastTest.txt