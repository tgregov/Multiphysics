#ifndef computeNorm_hpp_included
#define computeNorm_hpp_included

#include <Eigen/Dense>
#include "../mesh/Mesh.hpp"
#include "../params/Params.hpp"


double func(double x, double t, const std::vector<double>& coeffs, const std::vector<double>& fluxCoeffs);
/**
 * brief compute the analytical solution for a 1D gaussian wave, at position x and time t, for a Gaussian
 */

double func2(double x, double t, const std::vector<double>& coeffs, const std::vector<double>& fluxCoeffs);
/**
 * brief compute the analytical solution for a 1D gaussian wave, at position x and time t, for a parabola
 */


/**
 * \brief compute the L2 norm of the error
 * \param mesh The mesh representing the domain of interest
 * \param solverParams The structure in which the parameters of the solver are.
 * \param t time
 * \param u the solution computed by timeInteg
 * \return the value of the error
 */
void computeNorm(const Mesh& mesh, const SolverParams& solverParams, double t, Eigen::VectorXd u, double& errorL2, double& errorLinf);

#endif /* computeNorm_hpp */