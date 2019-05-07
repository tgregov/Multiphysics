#!/bin/bash
# Submission script for NIC4 
#SBATCH --job-name=Multiphysics
#SBATCH --time=10:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1024 # megabytes 
#SBATCH --partition=defq 
#SBATCH --output=outRK_part2.txt
##SBATCH --mail-user=joachim.marichal@student.uliege.be
#SBATCH --mail-type=ALL

module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=16

cd $HOME/Multiphysics
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_5.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2_5.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_4.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2_4.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_3.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2_3.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_2.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2_2.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_1.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2_1.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_05.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2_05.msh ./results/errorRK2.txt
srun ./build/bin/main ./Geometry/2D\ Rectangle/rectangle0_025.geo ./Geometry/2D\ Rectangle/newMesh.msh  ./Params/paramRK2.dat ./results/newResultsRK2_025.msh ./results/errorRK2.txt