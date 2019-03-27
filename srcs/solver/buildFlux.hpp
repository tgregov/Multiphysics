#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../mesh/Mesh2D.hpp"
#include "../params/ibvFunction.hpp"


/**
 * \brief Physical flux (allows to compute it at each node of the mesh).
 * \param fx Physical flux variable along x.
 * \param fy Physical flux variable along y.
 * \param C Parameter C of the Lax-Friedrichs numerical flux.
 * \param u Nodal vector of unknowns.
 */
void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, double& C,
			const Eigen::VectorXd& u);


/**
 * \brief Physical flux (allows to compute it at the boundary conditions).
 * \param fx Physical flux variable along x.
 * \param fy Physical flux variable along y.
 * \param u Nodal vector of unknowns.
 */
void flux(double& fx, double& fy, double u);


/**
 * \brief Function that allows to build the rhs of the DG method.
 * \param mesh The mesh of the problem.
 * \param I The rhs vector of the problem.
 * \param u Nodal vector of unknowns.
 * \param fx Physical flux variable along x, evaluated at each node.
 * \param fy Physical flux variable along y, evaluated at each node.
 * \param C Parameter C of the Lax-Friedrichs numerical flux.
 * \param factor Parameter that allows to get the weak or stong form.
 * \param numNodes Number of nodes within the mesh.
 * \param t Time.
 * \param boundaries Information about the boundaries.
 */
void buildFlux(const Mesh2D& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u,
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy, double C,
	double factor, unsigned int numNodes, double t,
	const std::map<std::string, bc>& boundaries);

#endif /* buildFlux_hpp */
