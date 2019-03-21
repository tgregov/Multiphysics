#include <iostream>
#include <string>
#include <Eigen/Sparse>
#include "./mesh/Mesh2D.hpp"
#include "./mesh/displayMesh.hpp"
#include "timeInteg.hpp"
#include "Solver.hpp"

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
    SolverParams solverParams;
    if(!loadSolverParams(std::string(argv[2]), solverParams))
        return -1;

    // load the mesh
    Mesh2D mesh;
    if(!readMesh2D(mesh, std::string(argv[1]), solverParams.spaceIntType, 
        solverParams.basisFuncType))
    {
        std::cerr   << "Something went wrong when reading mesh file: "
                    << argv[1] << std::endl;
        return -1;
    }
    // displayMesh(mesh);
    
    // time-integrate the numerical scheme
    if(!timeInteg(mesh, solverParams.timeIntType, solverParams.timeStep,
                  solverParams.nbrTimeSteps, solverParams.solverType, 
                  std::string(argv[1])))
    {
        std::cerr   << "Something went wrong when time integrating" << std::endl;
        return -1;
    }

    return 0;
}
