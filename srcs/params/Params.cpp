/**
 * \file Params.cpp
 * \brief Implementation of the required function to load a SolverParams struct from file.
 */

#include <iostream>
#include <fstream>
#include <cerrno> //Change with something more c++ ;-)
#include <cstring>
#include <vector>
#include "Params.hpp"
#include "../physics/boundaryConditions.hpp"
#include "../physics/fluxes.hpp"
#include "../physics/phiPsis.hpp"
#include "../physics/sources.hpp"

/**
 * \brief Wrapper of std::getline to have comments in file with //
 * \param file The file which has been previously opened.
 * \param line The computed line (without comments)
 */
static void getLine(std::ifstream& file, std::string& line)
{
    std::string tempLine;
    while(true) //be careful with this
    {
        std::getline(file, tempLine);
        if(tempLine.compare(0, 2, "//") != 0)
        {
            std::size_t pos;
            pos = tempLine.find("//");
            if(pos != std::string::npos)
            {
                while(true)
                {
                    if(tempLine.at(pos-1) == ' ' || tempLine.at(pos-1) == '\t')
                        tempLine.erase(pos-1, 1);
                    else
                    {
                        tempLine.erase(pos);
                        break;
                    }
                    pos = tempLine.find("//");
                }
            }

            line = tempLine;
            return;
        }
    }
}


/**
 * \brief Load boundary conditions from a file
 * \param paramFile The file which has been previously opened.
 * \param solverParams The structure in which the parameters are loaded.
 * \param fileName The name of the parameters file which has been opened (for debug output).
 * \return true if the loading succeeds, false otherwise.
 */
static bool handleBoundaryCondition(std::ifstream& paramFile, SolverParams& solverParams,
                                    const std::string& fileName)
{
    unsigned int nBC = 0;
    bool foundInitCond = false;

    while(true)
    {
        ibc tempCondition;
        std::string bcName;
        std::string bcType;
        std::string tempBcCoeff;
        std::vector<double> bcCoeff;

        if(paramFile.eof())
        {
            std::cout << "End of file reached." << std::endl;
            break;
        }
        getLine(paramFile, bcName);

        getLine(paramFile, bcType);
        if(bcType[0] != '\t') //Let's make them tabulated for easier reading
        {
            std::cerr << "Bad type format for BC  " << bcName
                      << " in parameter file " << fileName << std::endl;

            return false;
        }
        bcType.erase(0,1);

        bool error = false;
        if(solverParams.problemType == "transport")
        {
            if(bcType == "constant")
                tempCondition.ibcFunc = constant;

            else if(bcType == "sinusTransport")
                tempCondition.ibcFunc = sinusTransport;

            else if(bcType == "gaussianTransport")
                tempCondition.ibcFunc = gaussianTransport;

            else if(bcType == "freeTransport")
                tempCondition.ibcFunc = freeTransport;

            else if(bcType == "gaussian2DTransport")
                tempCondition.ibcFunc = gaussian2DTransport;

            else
                error = true;
        }
        else if(solverParams.problemType == "shallow")
        {
            if(bcType == "constant")
                tempCondition.ibcFunc = constant;

            else if(bcType == "affineShallow")
                tempCondition.ibcFunc = affineShallow;

            else if(bcType == "sinusShallow")
                tempCondition.ibcFunc = sinusShallow;

            else if(bcType == "sinusAffShallow")
                tempCondition.ibcFunc = sinusAffShallow;

            else if(bcType == "reflectShallow")
                tempCondition.ibcFunc = reflectShallow;

            else if(bcType == "gaussian2DShallow")
                tempCondition.ibcFunc = gaussian2DShallow;

            else if(bcType == "gaussian1DShallowX")
                tempCondition.ibcFunc = gaussian1DShallowX;

            else if(bcType == "gaussian1DShallowY")
                tempCondition.ibcFunc = gaussian1DShallowY;

            else if(bcType == "openShallow")
                tempCondition.ibcFunc = openShallow;

            else if(bcType == "openAffShallow")
                tempCondition.ibcFunc = openAffShallow;

            else
                error = true;

        }
        else if(solverParams.problemType == "shallowLin")
        {
            if(bcType == "constant")
                tempCondition.ibcFunc = constant;

            else if(bcType == "sinusShallowLin")
                tempCondition.ibcFunc = sinusShallowLin;

            else if(bcType == "reflectShallowLin")
                tempCondition.ibcFunc = reflectShallowLin;

            else if(bcType == "gaussian2DShallowLin")
                tempCondition.ibcFunc = gaussian2DShallowLin;

            else if(bcType == "gaussian1DShallowXLin")
                tempCondition.ibcFunc = gaussian1DShallowXLin;

            else if(bcType == "gaussian1DShallowYLin")
                tempCondition.ibcFunc = gaussian1DShallowYLin;

            else if(bcType == "openShallowLin")
                tempCondition.ibcFunc = openShallowLin;

            else
                error = true;
        }

        if(error)
        {
            std::cerr << "Unhandled boundary condition type " << bcType
                      << " for boundary " << std::endl << bcName
                      << " for problem type "
                      << solverParams.problemType << " in parameter file "
                      << fileName <<std::endl;

            return false;
        }

        getLine(paramFile, tempBcCoeff);
        if(tempBcCoeff[0] != '\t') //Let's make them tabulated for easier reading
        {
            std::cerr << "Bad coefficient format for BC  " << bcName
                      << " in parameter file " << fileName << std::endl;

            return false;
        }
        tempBcCoeff.erase(0,1);
        unsigned int precComaPos = -1;
        for(unsigned int i = 0 ; i < tempBcCoeff.size() ; ++i)
        {
            if(tempBcCoeff[i] == ',')
            {
                bcCoeff.push_back(std::stod(tempBcCoeff
                            .substr(precComaPos + 1, i - precComaPos - 1)));
                precComaPos = i;
            }
        }
        //At the end, still one push_back to do
        bcCoeff.push_back(std::stod(tempBcCoeff
                    .substr(precComaPos+1, tempBcCoeff.size() - precComaPos - 1)));

        tempCondition.coefficients = bcCoeff;

        if(bcName == "Init_Cond")
        {
            solverParams.initCondition = tempCondition;
            foundInitCond = true;
        }
        else
        {
            nBC++;
            solverParams.boundaryConditions[bcName] = tempCondition;
        }
    }

    if(!foundInitCond)
    {
        std::cerr << "No initial condition (Init_Cond field) found in parameter file"
                  << " " << fileName << std::endl;

        return false;
    }

    std::cout << "Initial condition present and " << nBC
              << " boundary conditions present in "<< std::endl
              << "file " << fileName << std::endl;

    return true;
}

//Documentation in .hpp
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
    getLine(paramFile, temp);

    if(temp.compare(0, 5, "Gauss") != 0
       || !(temp.substr(5, temp.size() - 5).find_first_not_of("0123456789")
            == std::string::npos))
    {
        std::cerr << "Unexpected space integration type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.spaceIntType = temp;

    temp.clear();
    getLine(paramFile, temp);

    if(!(temp == "Lagrange" || temp == "Isoparametric"))
    {
        std::cerr << "Unexpected basis function type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.basisFuncType = temp;

    temp.clear();
    getLine(paramFile, temp);

    if(!(temp == "RK1" || temp == "RK2" || temp == "RK3" || temp == "RK4"))
    {
        std::cerr << "Unexpected time integration type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.timeIntType = temp;

    temp.clear();
    getLine(paramFile, temp);

    if(!(temp == "strong" || temp == "weak"))
    {
        std::cerr << "Unexpected solver type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.solverType = temp;

    temp.clear();
    getLine(paramFile, temp);

    if(!(temp.find_first_not_of(".0123456789") == std::string::npos))
    {
        std::cerr << "Unexpected simulation time duration " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.simTime = std::stod(temp);

    temp.clear();
    getLine(paramFile, temp);

    if(!(temp.find_first_not_of(".0123456789") == std::string::npos))
    {
        std::cerr << "Unexpected time step " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    solverParams.timeStep = std::stod(temp);

    temp.clear();
    getLine(paramFile, temp);

    if(!(temp.find_first_not_of(".0123456789") == std::string::npos))
    {
        std::cerr << "Unexpected time between data writing " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }
    solverParams.simTimeDtWrite = std::stod(temp);

    temp.clear();
    getLine(paramFile, temp);
    if(temp == "shallow")
    {
        solverParams.problemType = temp;
        solverParams.nUnknowns = 3;
        solverParams.flux = fluxShallow;
    }
    else if(temp == "transport")
    {
        solverParams.problemType = temp;
        solverParams.nUnknowns = 1;
        solverParams.flux = fluxTransport;
    }
    else if(temp == "shallowLin")
    {
        solverParams.problemType = temp;
        solverParams.nUnknowns = 3;
        solverParams.flux = fluxShallowLin;
    }
    else
    {
        std::cerr << "Unexpected problem type " << temp
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }

    temp.clear();
    getLine(paramFile, temp);
    unsigned int precComaPos = -1;
    for(unsigned int i = 0 ; i < temp.size() ; ++i)
    {
        if(temp[i] == ',')
        {
            solverParams.whatToWrite.push_back(
                (temp.substr(precComaPos + 1, i - precComaPos - 1))
                 == "1" ? true : false);
            precComaPos = i;
        }
    }
    //At the end, still one push_back to do
    solverParams.whatToWrite.push_back(
        (temp.substr(precComaPos+1, temp.size() - precComaPos - 1))
        == "1" ? true : false);

    bool error = false;
    if(solverParams.problemType == "shallow"
        || solverParams.problemType == "shallowLin")
    {
        if(solverParams.whatToWrite.size() !=5)
            error = true;
        else
        {
            solverParams.write = writeShallow;
        }
    }
    else if(solverParams.problemType == "transport")
    {
        if(solverParams.whatToWrite.size() !=1)
            error = true;
        else
        {
            solverParams.write = writeTransport;
        }
    }
    if(error)
    {
        std::cerr << "Unexpected number of writing parameters ("
                  << solverParams.whatToWrite.size() << ") for problem type "
                  << solverParams.problemType << std::endl;

        paramFile.close();

        return false;
    }

    solverParams.viewTags.resize(solverParams.whatToWrite.size());

    temp.clear();
    getLine(paramFile, temp);
    error = false;
    if(solverParams.problemType == "shallow")
    {
        if(temp == "LF")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = LFShallow;
        }
        else if(temp == "Roe")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = Roe;
        }
        else if(temp == "mean")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = mean;
        }
        else
            error = true;
    }
    else if(solverParams.problemType == "shallowLin")
    {
        if(temp == "LF")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = LFShallowLin;
        }
        else if(temp == "Roe")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = RoeLin;
        }
        else if(temp == "mean")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = mean;
        }
        else
            error = true;
    }
    else if(solverParams.problemType == "transport")
    {
        if(temp == "LF")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = LFTransport;
        }
        else if(temp == "mean")
        {
            solverParams.fluxType = temp;
            solverParams.phiPsi = mean;
        }
        else
            error = true;
    }
    if(error)
    {
        std::cerr << "Unexpected flux type type " << temp
                  << " for problem type " << solverParams.problemType
                  << " in parameter file " << fileName << std::endl;

        paramFile.close();
        return false;
    }


    temp.clear();
    getLine(paramFile, temp);
    precComaPos = -1;
    for(unsigned int i = 0 ; i < temp.size() ; ++i)
    {
        if(temp[i] == ',')
        {
            solverParams.fluxCoeffs.push_back(
                    std::stod(temp.substr(precComaPos + 1, i - precComaPos - 1)));
            precComaPos = i;
        }
    }
    //At the end, still one push_back to do
    solverParams.fluxCoeffs.push_back(
            std::stod(temp.substr(precComaPos+1, temp.size() - precComaPos - 1)));

    error = false;
    if(solverParams.problemType == "shallow")
    {
        if(solverParams.fluxCoeffs.size() !=1)
            error = true;
    }
    else if(solverParams.problemType == "shallowLin")
    {
        if(solverParams.fluxCoeffs.size() != 2)
            error = true;

    }
    if(solverParams.problemType == "transport")
    {
        if(solverParams.fluxCoeffs.size() != 2)
            error = true;

    }
    if(error)
    {
        std::cerr << "Unexpected number of flux coefficients ("
                  << solverParams.fluxCoeffs.size() << ") for problem type "
                  << solverParams.problemType << std::endl;

        paramFile.close();

        return false;
    }

    temp.clear();
    getLine(paramFile, temp);
    if(temp == "no")
    {
        solverParams.IsSourceTerms = false;
        solverParams.sourceType = temp;
    }
    else if(temp == "sourceShallowCstGradCstFrict"
            && solverParams.problemType == "shallow")
    {
        solverParams.IsSourceTerms = true;
        solverParams.sourceType = temp;
        solverParams.sourceTerm = sourceShallowCstGradCstFrict;
    }
    else if(temp == "sourceShallowCstGradQuadFrict"
            && solverParams.problemType == "shallow")
    {
        solverParams.IsSourceTerms = true;
        solverParams.sourceType = temp;
        solverParams.sourceTerm = sourceShallowCstGradQuadFrict;
    }
    else if(temp == "shallowLinCst" && solverParams.problemType == "shallowLin")
    {
        solverParams.IsSourceTerms = true;
        solverParams.sourceType = temp;
        solverParams.sourceTerm = sourceShallowLinCst;
    }
    else
    {
        std::cerr << "Unexpected source function ("
                  << solverParams.sourceType << ") for problem type "
                  << solverParams.problemType << std::endl;

        paramFile.close();
        return false;
    }

    temp.clear();
    getLine(paramFile, temp);
    precComaPos = -1;
    for(unsigned int i = 0 ; i < temp.size() ; ++i)
    {
        if(temp[i] == ',')
        {
            solverParams.sourceCoeffs.push_back(
                std::stod(temp.substr(precComaPos + 1, i - precComaPos - 1)));
            precComaPos = i;
        }
    }
    //At the end, still one push_back to do
    solverParams.sourceCoeffs.push_back(
        std::stod(temp.substr(precComaPos+1, temp.size() - precComaPos - 1)));


    error = false;
    if(solverParams.sourceType == "sourceShallowCstGradCstFrict")
    {
        if(solverParams.IsSourceTerms)
        {
            if(solverParams.sourceCoeffs.size() != 5)
                error = true;
        }
    }
    if(solverParams.sourceType == "sourceShallowCstGradQuadFrict")
    {
        if(solverParams.IsSourceTerms)
        {
            if(solverParams.sourceCoeffs.size() != 5)
                error = true;
        }
    }
    else if(solverParams.sourceType == "shallowLinCst")
    {
        if(solverParams.IsSourceTerms)
        {
            if(solverParams.sourceCoeffs.size() != 1)
                error = true;
        }
    }

    if(error)
    {
        std::cerr << "Unexpected number of source terms coefficients ("
                      << solverParams.sourceCoeffs.size() << ") for problem type "
                      << solverParams.problemType << std::endl;

        paramFile.close();

        return false;
    }

    if(!handleBoundaryCondition(paramFile, solverParams, fileName))
    {
        paramFile.close();
        return false;
    }

    paramFile.close();

    // display the parameters
    std::cout   << "Number of Gauss points: " << solverParams.spaceIntType
                << std::endl
                << "Type of basis function: " << solverParams.basisFuncType
                << std::endl
                << "Time intgeration scheme: " << solverParams.timeIntType
                << std::endl
                << "Formulation type: " << solverParams.solverType
                << std::endl
                << "Simulation time duration: " << solverParams.simTime << "s"
                << std::endl
                << "Time step: " << solverParams.timeStep << "s"
                << std::endl
                << "Time between data writing: " << solverParams.simTimeDtWrite
                << "s"
                << std::endl
                << "Problem type: " << solverParams.problemType
                << std::endl
                << "Source terms: " << solverParams.sourceType
                << std::endl
                << "Numerical Flux: " << solverParams.fluxType
                << std::endl;

    return true;
}
