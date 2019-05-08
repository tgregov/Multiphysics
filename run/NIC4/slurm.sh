#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=01:00:00 # hh:mm:ss
#
#SBATCH --ntasks=4 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=out.txt

module load gcc/4.9.2
module load openmpi/1.8.4/gcc-4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=16
export OMP_CANCELLATION=true
export MPI_NUM_THREADS=4

cd $HOME/Multiphysics
mpirun --bind-to none  -np $MPI_NUM_THREADS ./build/bin/main ./Geometry/2D\ Rectangle/rectangle.msh ./Params/param.dat
