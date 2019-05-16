#!/bin/bash
# Submission script for VEGA 
#SBATCH --job-name=Antarctic
#SBATCH --time=15:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=4000 # megabytes 
#SBATCH --partition=defq 
#
#SBATCH --comment=Multiphysics 
#SBATCH --output=outAnt.txt
#
#SBATCH --mail-user=joachim.marichal@student.uliege.be
#SBATCH --mail-type=ALL

module load gcc/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_NUM_THREADS=64
export OMP_CANCELLATION=true

cd $HOME/Multiphysics
srun ./build/bin/main ./Geometry/mer/antarctique3.geo ./Geometry/mer/antarctique3Phys.msh ./Params/paramIceM.dat ./results/resultsAntarctic.msh


