#include <iostream>
#include <string>
#if defined(_OPENMP)
    #include <cstdlib>
    #include <omp.h>
#endif
#include <Eigen/Core>
#include "mesh/Mesh.hpp"
#include "mesh/displayMesh.hpp"
#include "solver/timeInteg.hpp"
#include "params/Params.hpp"
#include "error/computeError.hpp"

int main(int argc, char **argv)
{

    // check that the file format is valid
    if (argc < 3)
    {
        std::cerr   << "Usage: " << argv[0] << " file.msh " << " param.dat "
                    <<  std::endl;
        return 1;
    }

    // load the solver parameters
    std::cout   << "================================================================"
                << std::endl
                << "                   LOADING THE SOLVER PARAMETERS                "
                << std::endl
                << "================================================================"
                << std::endl;

    SolverParams solverParams;
    if(!loadSolverParams(std::string(argv[2]), solverParams))
        return -1;

    #if defined(_OPENMP)
        unsigned int n = std::atoi(std::getenv("OMP_NUM_THREADS"));
        omp_set_num_threads(n);
        Eigen::setNbThreads(1);
        std::cout << "Number of threads: " << n << std::endl;;
    #endif

    // load the mesh
    std::cout   << "================================================================"
                << std::endl
                << "                       LOADING THE MESH                         "
                << std::endl
                << "================================================================"
                << std::endl;
    Mesh mesh;
    if(!readMesh(mesh, std::string(argv[1]), solverParams.spaceIntType,
        solverParams.basisFuncType))
    {
        std::cerr   << "Something went wrong when reading mesh file: "
                    << argv[1] << std::endl;
        return -1;
    }

   // displayMesh(mesh);
   std::cout   << "================================================================"
                << std::endl
                << "                     EXECUTING THE SOLVER                       "
                << std::endl
                << "================================================================"
                << std::endl;
    if(!timeInteg(mesh, solverParams, std::string(argv[1]), std::string(argv[3])))
    {
        std::cerr   << "Something went wrong when time integrating" << std::endl;
        return -1;
    }

      std::cout   << "================================================================"
                << std::endl
                << "                       ERROR COMPUTATION                           "
                << std::endl
                << "==================================================================="
                << std::endl;
    if(!computeError(mesh, std::string(argv[1]), std::string(argv[3])))
    {
        std::cerr << "Something went wrong when computing the error" << std::endl;
        return -1;
    }

    return 0;
}