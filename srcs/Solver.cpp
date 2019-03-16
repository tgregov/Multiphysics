#include <iostream>
#include <fstream>
#include <cerrno> //Change with something more c++ ;-)
#include <cstring>
#include "Solver.hpp"

bool loadSolverParams(const std::string& fileName, SolverParams& solverParams)
{
    std::ifstream paramFile(fileName);

    if(!paramFile.is_open())
    {
        std::cerr << "Something went wrong when trying to read the file "
                  << fileName << std::endl
                  << std::strerror(errno) << std::endl;

        paramFile.close();
        return false;
    }

    std::string temp;
    std::getline(paramFile, temp);

    if(temp.find("Gauss") != 0)
    {
        std::cerr << "Unexpected space integration type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.spaceIntType = temp;

    temp.clear();
    std::getline(paramFile, temp);

    if(!(temp == "Lagrange" || temp == "Isoparametric"))
    {
        std::cerr << "Unexpected basis function type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.basisFuncType = temp;

    temp.clear();
    std::getline(paramFile, temp);

    if(!(temp == "RK1" || temp == "RK4"))
    {
        std::cerr << "Unexpected time integration type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.timeIntType = temp;

    temp.clear();
    std::getline(paramFile, temp);

    if(!(temp == "strong" || temp == "weak"))
    {
        std::cerr << "Unexpected solver type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.solverType = temp;

    temp.clear();
    std::getline(paramFile, temp);

    if(!(temp.find_first_not_of("0123456789") == std::string::npos)) //To improve
    {
        std::cerr << "Unexpected number of time steps " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.nbrTimeSteps = std::stoi(temp);

    temp.clear();
    std::getline(paramFile, temp);

    if(!(temp.find_first_not_of(".0123456789") == std::string::npos)) //To improve
    {
        std::cerr << "Unexpected time step " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.timeStep = std::stod(temp);

    paramFile.close();
    return true;
}
