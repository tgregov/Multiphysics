#include <iostream>
#include <string>
#if defined(_OPENMP)
    #include <cstdlib>
    #include <omp.h>
#endif
#include <Eigen/Core>
#include <mpi.h>
#include "mesh/Mesh.hpp"
#include "mesh/displayMesh.hpp"
#include "solver/timeInteg.hpp"
#include "params/Params.hpp"

int main(int argc, char **argv)
{

    // check that the file format is valid
    int numberProc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberProc);

    if (argc < 3)
    {
        std::cerr   << "Usage: " << argv[0] << " file.msh " << " param.dat "
                    <<  std::endl;
        return 1;
    }

    // load the solver parameters
    if(rank == 0)
    {
        std::cout   << "================================================================"
                    << std::endl
                    << "                   LOADING THE SOLVER PARAMETERS                "
                    << std::endl
                    << "================================================================"
                    << std::endl;
    }


    SolverParams solverParams;
    if(!loadSolverParams(std::string(argv[2]), solverParams, rank))
    {
        MPI_Finalize();
        return -1;
    }

    #if defined(_OPENMP)
        unsigned int n = std::atoi(std::getenv("OMP_NUM_THREADS"));
        omp_set_num_threads(n);
        Eigen::setNbThreads(1);
        if(rank == 0)
            std::cout << "Number of threads: " << n << std::endl;
    #endif

    if(rank == 0)
    {
        std::cout << "Number of MPI threads: " << numberProc << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // load the mesh
    if(rank == 0)
    {
        std::cout   << "================================================================"
                    << std::endl
                    << "                       LOADING THE MESH                         "
                    << std::endl
                    << "================================================================"
                    << std::endl;
    }

    Mesh mesh;

    auto startTime = MPI_Wtime();
    if(!readMesh(mesh, std::string(argv[1]), solverParams.spaceIntType,
        solverParams.basisFuncType, rank))
    {
        std::cerr   << "Something went wrong when reading mesh file: "
                    << argv[1] << std::endl;

        MPI_Finalize();
        return -1;
    }
    auto endTime =  MPI_Wtime();
    auto ellapsedTime = endTime - startTime;

    if(rank == 0)
    {
        std::cout << "Ellapsed time for mesh reading: "
                  << ellapsedTime
                  << " s" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // displayMesh(mesh);
    if(rank == 0)
    {
        std::cout   << "================================================================"
                    << std::endl
                    << "                     EXECUTING THE SOLVER                       "
                    << std::endl
                    << "================================================================"
                    << std::endl;
    }

    startTime =  MPI_Wtime();
    if(!timeInteg(mesh, solverParams, std::string(argv[1]), rank, numberProc))
    {
        std::cerr   << "Something went wrong when time integrating" << std::endl;
        MPI_Finalize();
        return -1;
    }
    endTime =  MPI_Wtime();
    ellapsedTime = endTime - startTime;

    if(rank == 0)
    {
        std::cout << "Ellapsed time for time integration: "
                  << ellapsedTime
                  << " s" << std::endl;
    }

    MPI_Finalize();

    return 0;
}
