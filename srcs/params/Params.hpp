#ifndef Params_hpp_included
#define Params_hpp_included

#include <vector>
#include <string>
#include <functional>
#include "../solver/field.hpp"
#include "../mesh/Mesh.hpp"
#include "../physics/writers.hpp"
#include "../physics/ibcFunction.hpp"
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

    std::map<std::string, ibc> boundaryConditions;  /**< Map of the problem's boundary condition*/
    ibc initCondition;                              /**< Initial condition*/

    std::string problemType;     /**< Equations to solve (transport, shallow, ...)*/

    std::function<void(Field& field,
                       PartialField& partialField,
                       const SolverParams& solverParams,
                       bool boundary)> flux;    /**< Pointer to the flux function
                                            (represents the physics of the problem)*/

    std::vector<double> fluxCoeffs; /**< Coefficient of the physical flux*/

    unsigned short nUnknowns;    /**< Number of scalar unknowns
                                      (determined by the type of problem)*/

    std::string fluxType;        /**< Type of numerical flux
                                      (mean, Lax-Friedirichs, Roe, ...)*/

    std::function<void(const Edge& edge, Field& field, PartialField& partialField,
                       unsigned int j, double factor, bool boundary,
                       unsigned int indexJ, unsigned int indexFrontJ,
                       const SolverParams& solverParams)> phiPsi; /**< Pointer to the
                       rhs function (phi or psi depending of the type of scheme)*/

    bool IsSourceTerms;                 /**< (De)activate source terms computation*/
    std::string sourceType;             /**< Denotes the type of source terms*/
    std::vector<double> sourceCoeffs;   /**< Coefficient of the source terms*/
    std::function<void(Field& field, const SolverParams& solverParams)> sourceTerm;/**< Pointer to the source terms function*/


    std::vector<bool> whatToWrite; /**< Vector of boolean denoting
                                        what will be written (problem dependent)*/
    std::vector<int> viewTags;  /**< Store the view tag of what will be written*/

    std::function<void(std::vector<std::vector<double>>& uDisplay,
                  const std::vector<unsigned int>& elementNumNodes,
                  const std::vector<int>& elementTags, const std::string& modelName,
                  unsigned int nbreStep, double t, const Field& field,
                  const std::vector<double>& fluxCoeffs,
                  const std::vector<bool>& whatToWrite,
                  std::vector<int>& viewTags)> write; /**< Pointer to the function
                                                           which will write*/
};

/**
 * \brief Load solver parameters from file in a structure
 * \param fileName The name of the parameters file to load.
 * \param solverParams The structure in which the parameters are loaded.
 * \return true if the loading succeeds, false otherwise.
 */
bool loadSolverParams(const std::string& fileName, SolverParams& solverParams);

#endif // Params_hpp_included
