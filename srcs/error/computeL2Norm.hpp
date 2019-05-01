#ifndef computeL2Norm_hpp_included
#define computeL2Norm_hpp_included

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
double computeL2Norm(const Mesh& mesh, const SolverParams& solverParams, double t, Eigen::VectorXd u);

#endif /* computeL2Norm_hpp */
