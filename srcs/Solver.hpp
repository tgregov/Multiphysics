#ifndef Solver_hpp_included
#define Solver_hpp_included

#include <string>
#include "bcFunction.hpp"

/**
 * \struct SolverParams
 * \brief Parameters of the method.
 */
struct SolverParams
{
    std::string spaceIntType;   /**< Number of points for Gauss integration
                                   (format Gaussx, x the number of points) */
    std::string basisFuncType;  /**< Type of basis functions (Lagrange or Isoparametric */
    std::string timeIntType;    /**< Runge-Kutta time integration type (RK1 or RK4) */
    std::string solverType;     /**< Solver form (strong or weak) */

    unsigned int simTime;       /**< Simulation time duration */
    double timeStep;            /**< Time steps for the simulation */

    std::map<std::string, bc> boundaryConditions;
};

/**
 * \brief Load solver parameters from file in a structure
 * \param fileName The name of the parameters file to load.
 * \param solverParams The structure in which the parameters are loaded.
 * \return true if the loading succeeds, false otherwise.
 */
bool loadSolverParams(const std::string& fileName, SolverParams& solverParams);

#endif // Solver_hpp_included
