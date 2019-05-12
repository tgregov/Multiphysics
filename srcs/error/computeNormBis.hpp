#ifndef computeNormBis_hpp_included
#define computeNormBis_hpp_included

#include <Eigen/Dense>
#include "../mesh/Mesh.hpp"
#include "../params/Params.hpp"

/**
 * \brief compute the L2 norm of the error
 * \param mesh The mesh representing the domain of interest
 * \param solverParams The structure in which the parameters of the solver are.
 * \param t time
 * \param u the solution computed by timeInteg
 * \return the value of the error
 */
void computeNormBis(const Mesh& mesh, const SolverParams& solverParams, double t, Eigen::VectorXd u, double& errorL2, double& errorLinf);

#endif /* computeNormBis_hpp */