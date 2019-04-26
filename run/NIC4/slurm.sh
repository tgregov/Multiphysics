#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=01:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=out.txt

module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=10

cd $HOME/Multiphysics
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle.msh ./Params/param.dat