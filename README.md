# MATH 0471: Multiphysics integrated computational project 
## Goal
Development of the [Discontinuous Galerkin method](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) (DG) applied to 2D [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations).  
[Link](http://www.montefiore.ulg.ac.be/~geuzaine/MATH0471/enonce2019.pdf) of the problem statement.

## Useful reminders
[Link](https://github.com/tgregov/Multiphysics/wiki) to the wiki.

## Build status 
Current status: [![Build Status](https://travis-ci.org/tgregov/Multiphysics.svg?branch=master)](https://travis-ci.org/tgregov/Multiphysics)

## Assumptions
### On the mesh:
* The elements edges are straight (constant normal vector over it)

## Compilation procedure
### On NIC4
Connect to NIC4 (using SSH for instance). Then, clone the repositoy: 
```bash
git clone https://github.com/tgregov/Multiphysics
```
Move to the code repository:
```bash
cd ./Multiphysics
```
Automatically build the code:
```bash
. ./build_nic4.sh 
```
Submit a batch file:
```bash
cd ./run/NIC4
sbatch slurm.sh
```