# MATH 0471: Multiphysics integrated computational project 
 [![Build Status](https://travis-ci.org/tgregov/Multiphysics.svg?branch=master)](https://travis-ci.org/tgregov/Multiphysics) [![Maintenance](https://img.shields.io/badge/Version-1.1.1-e67e22.svg)](https://github.com/tgregov/Multiphysics/releases/tag/1.1.1) 
 
## Description
This repository contains the implementation of the [Discontinuous Galerkin method](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) (DG) applied to systems of possibly coupled non-linear hyperbolic equations, i.e. to physical conservative systems. The geometry is currently assumed to be two-dimensional. The targeted physical application was chosen to be the [shallow water equations](https://en.wikipedia.org/wiki/Shallow_water_equations), but any new system of equations can easily be implemented, by adding new physical fluxes or sources.
 
The time-integration process is explicit, and uses Runge-Kutta algorithms, up to the fourth order. The mesh is handled using [Gmsh](http://gmsh.info/), and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is used for the linear algebra operations. Please note that edges of the elements in the mesh are assumed to be straight.

The code is currently parallized locally using [OpenMP](https://www.openmp.org/), such that it can be run on HPC clusters.
 
This project was part of the 2018-2019 *Multiphysics integrated computational project* course. Here is a [link](http://www.montefiore.ulg.ac.be/~geuzaine/MATH0471/enonce2019.pdf) of the problem statement. The authors are [*ImperatorS79*](https://github.com/ImperatorS79), [*tgregov*](https://github.com/tgregov/) and [*jmarichal*](https://github.com/jmarichal).

## Useful information
Some useful information can be found on the [wiki](https://github.com/tgregov/Multiphysics/wiki). In particular, a small (uncomplete) description of the physical parameters that can be introduced in the context of the shallow water equations is [provided](https://github.com/tgregov/Multiphysics/wiki/How-to-change-the-simulation-parameters-%3F).

## Compilation procedure
The project use [CMake](https://cmake.org/) to manage the build system. The `run` directory contains some predefined configuration for various platforms.

**Do not forget to install the [Gmsh 4.4.0 SDK](http://gmsh.info/) as well as the [Eigen 3.3.7 SDK](http://eigen.tuxfamily.org/index.php?title=Main_Page) on your system.**

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
