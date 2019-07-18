# MATH 0471: Multiphysics integrated computational project 
 [![Build Status](https://travis-ci.org/tgregov/Multiphysics.svg?branch=master)](https://travis-ci.org/tgregov/Multiphysics) [![Maintenance](https://img.shields.io/badge/Version-1.1.1-e67e22.svg)](https://github.com/tgregov/Multiphysics/releases/tag/1.1.1) 
 
Implementation of the [Discontinuous Galerkin method](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) (DG) applied to 2D [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations) using [gmsh](http://gmsh.info/) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page).  

[Link](http://www.montefiore.ulg.ac.be/~geuzaine/MATH0471/enonce2019.pdf) of the problem statement.

## Useful reminders
[Link](https://github.com/tgregov/Multiphysics/wiki) to the wiki.

## Assumptions
### On the mesh:
* The elements edges are straight (constant normal vector over it)

## Compilation procedure
The project use system to manage the build system. The `run` directory contains already some predefined configuration for various platform.

Do not forget to install gmsh 4.4.0 SDK and Eigen 3.3.7 SDK on your system.

### On NIC4/VEGA
Connect to NIC4/VEGA (using SSH for instance). Then, clone the repositoy: 
```bash
git clone https://github.com/tgregov/Multiphysics
```
Move to the code repository:
```bash
cd ./Multiphysics
```
Automatically build the code (`<cluster>` is either `vega` or `nic4`):
```bash
. ./build_cluster.sh <cluster>
```
Submit a batch file for one of the simulations
```bash
cd ../simulations
sbatch runObstacleSquare.sh
```
