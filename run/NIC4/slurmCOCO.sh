#!/bin/bash
# Submission script for VEGA
#SBATCH --job-name=Shallow
#SBATCH --time=16:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2048 # megabytes 
#SBATCH --partition=defq 
#
#SBATCH --comment=Multiphysics 
#SBATCH --output=outWeakScalling.txt

module load GCC/4.9.2
export CC=gcc
export CXX=g++
export FC=gfortran
export OMP_CANCELLATION=true

cd $HOME/Multiphysics

export OMP_NUM_THREADS=1
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle1.msh ./run/NIC4/weak/param1.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=2
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle2.msh ./run/NIC4/weak/param2.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=3
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle3.msh ./run/NIC4/weak/param3.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=4
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle4.msh ./run/NIC4/weak/param4.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=5
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle5.msh ./run/NIC4/weak/param5.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=6
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle6.msh ./run/NIC4/weak/param6.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=7
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle7.msh ./run/NIC4/weak/param7.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=8
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle8.msh ./run/NIC4/weak/param8.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=9
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle9.msh ./run/NIC4/weak/param9.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=10
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle10.msh ./run/NIC4/weak/param10.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=11
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle11.msh ./run/NIC4/weak/param11.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=12
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle12.msh ./run/NIC4/weak/param12.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=13
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle13.msh ./run/NIC4/weak/param13.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=14
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle14.msh ./run/NIC4/weak/param14.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=15
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle15.msh ./run/NIC4/weak/param15.dat ./run/NIC4/weak/result.msh

export OMP_NUM_THREADS=16
srun ./build/bin/main ./run/NIC4/weak/rectangle.geo ./run/NIC4/weak/rectangle16.msh ./run/NIC4/weak/param16.dat ./run/NIC4/weak/result.msh