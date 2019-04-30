#ifndef computeL2Norm_hpp_included
#define computeL2Norm_hpp_included

#include <Eigen/Dense>
#include "../mesh/Mesh.hpp"


/**
 * \brief Time integrate the equations (DG-FEM)
 * \param mesh The mesh representing the domain of interest
 * \param solverParams The structure in which the parameters of the solver are.
 * \param fileName The name of the file containing the mesh.
 * \return true if time integration happened without problems, false otherwise.
 */
double computeL2Norm(const Mesh& mesh, double t, Eigen::VectorXd u);

#endif /* computeError2_hpp */
