#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=10:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=outLF.txt
##SBATCH --mail-user=joachim.marichal@student.uliege.be
#SBATCH --mail-type=ALL

module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=16

cd $HOME/Multiphysics
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_5.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramLF.dat ./results/newResultsLF5.msh ./results/errorLF.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_4.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramLF.dat ./results/newResultsLF4.msh ./results/errorLF.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_3.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramLF.dat ./results/newResultsLF3.msh ./results/errorLF.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_2.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramLF.dat ./results/newResultsLF2.msh ./results/errorLF.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_1.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramLF.dat ./results/newResultsLF1.msh ./results/errorLF.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_05.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramLF.dat ./results/newResultsLF05.msh ./results/errorLF.txt