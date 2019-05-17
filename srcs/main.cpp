#include <chrono>
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

/**
 * @param  argv[1] .msh file that contains the mesh.
 * @param  argv[2] .dat file that contains the parameters.
 * @param  argv[3] name of the .msh file that will contain the results. 
 */
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

    // set the desired number of OpenMP threads
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
    auto startTime = std::chrono::high_resolution_clock::now();
    if(!readMesh(mesh, std::string(argv[1]), solverParams.spaceIntType,
        solverParams.basisFuncType))
    {
        std::cerr   << "Something went wrong when reading mesh file: "
                    << argv[1] << std::endl;
        return -1;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    auto ellapsedTime = 
        std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "Ellapsed time for mesh reading: "
              << static_cast<double>(ellapsedTime.count())/1000.0
              << " s" << std::endl;

   // displayMesh(mesh);
    std::cout   << "================================================================"
                << std::endl
                << "                     EXECUTING THE SOLVER                       "
                << std::endl
                << "================================================================"
                << std::endl;

    startTime = std::chrono::high_resolution_clock::now();
    if(!timeInteg(mesh, solverParams, std::string(argv[1]),std::string(argv[3])))
    {
        std::cerr   << "Something went wrong when time integrating" << std::endl;
        return -1;
    }
    endTime = std::chrono::high_resolution_clock::now();
    ellapsedTime =
        std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);

    std::cout << "Ellapsed time for time integration: "
              << static_cast<double>(ellapsedTime.count())/1000.0
              << " s" << std::endl;

    return 0;
}
