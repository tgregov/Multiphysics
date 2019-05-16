#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=01:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=outRK2.txt
##SBATCH --mail-user=joachim.marichal@student.uliege.be
#SBATCH --mail-type=ALL

module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=16

cd $HOME/Multiphysics
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangleLong0_5.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangleLong0_4.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangleLong0_3.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangleLong0_2.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangleLong0_1.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangleLong0_05.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangleLong0_025.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2.msh ./results/errorRK2.txt