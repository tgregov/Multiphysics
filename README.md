# MATH 0471: Multiphysics integrated computational project 
## Goal
Development of the Discontinuous Galerkin method (DG) applied to a particular physical situation (to be determined).  
[Link](http://www.montefiore.ulg.ac.be/~geuzaine/MATH0471/enonce2019.pdf) of the problem statement.

## Useful reminders
[Link](https://github.com/tgregov/Multiphysics/wiki) to the wiki.

## Build status 
Current status: [![Build Status](https://travis-ci.org/tgregov/Multiphysics.svg?branch=master)](https://travis-ci.org/tgregov/Multiphysics)

## Assumptions
### On the field:
* Scalar field

### On the mesh:
* 2D mesh
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
Run the code:
```bash
cd ..
./build/bin/main ./Geometry/2D\ Rectangle/rectangle.msh ./Params/param.dat
```


