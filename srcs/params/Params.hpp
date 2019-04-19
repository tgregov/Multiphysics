#ifndef Params_hpp_included
#define Params_hpp_included

#include <vector>
#include <string>
#include <functional>
#include "ibvFunction.hpp"
#include "../solver/field.hpp"
#include "../mesh/Mesh.hpp"

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

    double simTime;             /**< Simulation time duration */
    double timeStep;            /**< Time steps for the simulation */
    double simTimeDtWrite;      /**< Time between two data writings */

    std::map<std::string, ibc> boundaryConditions;
    ibc initCondition;

    std::string problemType;     /**< Equations to solve (transport, shallow, ...)*/

    std::function<void(Field& field,
                       const SolverParams& solverParams,
                       bool boundary)> flux;

    std::vector<double> fluxCoeffs; /**< Coefficient of the physical flux*/

    unsigned short nUnknowns;    /**< Number of scalar unknowns (determined by the type of problem)*/

    std::string fluxType;        /**< Type of numerical flux (mean, Lax-Friedirichs, Roe, ...)*/

    std::function<void(const Edge& edge, Field& field, unsigned int j,
                       double factor, bool boundary, unsigned int indexJ,
                       unsigned int indexFrontJ, const SolverParams& solverParams)> phiPsy;

};

/**
 * \brief Load solver parameters from file in a structure
 * \param fileName The name of the parameters file to load.
 * \param solverParams The structure in which the parameters are loaded.
 * \return true if the loading succeeds, false otherwise.
 */
bool loadSolverParams(const std::string& fileName, SolverParams& solverParams);

#endif // Params_hpp_included
