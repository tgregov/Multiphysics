#ifndef SOLVER_HPP_INCLUDED
#define SOLVER_HPP_INCLUDED

#include <string>

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

    unsigned int nbrTimeSteps;  /**< Number of time steps for the simulation */  //Change it for duration ?
    double timeStep;            /**< Time steps for the simulation */
};

/**
 * \brief Load solver parameters from file in a structure
 * \param fileName The name of the parameters file to load.
 * \param solverParams The structure in which the parameters are loaded.
 */
bool loadSolverParams(const std::string& fileName, SolverParams& solverParams);

#endif // SOLVER_HPP_INCLUDED
