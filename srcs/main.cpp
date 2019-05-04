#include <iostream>
#include <fstream>
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
#include "generateMesh.hpp"

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
    if(!loadSolverParams(std::string(argv[3]), solverParams))
        return -1;

    #if defined(_OPENMP)
        unsigned int n = std::atoi(std::getenv("OMP_NUM_THREADS"));
        omp_set_num_threads(n);
        Eigen::setNbThreads(1);
        std::cout << "Number of threads: " << n << std::endl;;
    #endif

    std::ofstream file1("errorVSflux.txt", std::ios::out | std::ios::app); 
    if (file1)
    {
            file1 << "mesh:" << "\t" << argv[1] << std::endl;
            file1 << "param:" << "\t" << argv[3] << std::endl;
            file1 << "ORDER:" << "\t" << "ERROR:" << std::endl;

    }else
            std::cerr << "Impossible d'ouvrir le fichier !" << std::endl;

    unsigned int order;
    double errorValue = 0;
    for (order = 1 ; order <= 5 ; order ++)
    {
        generateMesh(argv[1], argv[2], order);  

        std::cout   << "================================================================"
                    << std::endl
                    << "                       ORDER:                                   "
                    << std::endl
                    << "                            "
                    << order
                    << std::endl
                    << "================================================================"
                    << std::endl; 
        // load the mesh
        std::cout   << "================================================================"
                    << std::endl
                    << "                       LOADING THE MESH                         "
                    << std::endl
                    << "================================================================"
                    << std::endl;
        Mesh mesh;
        if(!readMesh(mesh, std::string(argv[2]), solverParams.spaceIntType,
            solverParams.basisFuncType))
        {
            std::cerr   << "Something went wrong when reading mesh file: "
                        << argv[2] << std::endl;
            return -1;
        }

       // displayMesh(mesh);
       std::cout   << "================================================================"
                    << std::endl
                    << "                     EXECUTING THE SOLVER                       "
                    << std::endl
                    << "================================================================"
                    << std::endl;
        if(!timeInteg(mesh, solverParams, std::string(argv[2]), std::string(argv[4])))
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
        if(!computeError(mesh, solverParams, std::string(argv[2]), std::string(argv[4]), errorValue))
        {
            std::cerr << "Something went wrong when computing the error" << std::endl;
            return -1;
        }
        file1 << order << "\t" << errorValue << std::endl;
    }


    file1.close();
    return 0;
}