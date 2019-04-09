#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../mesh/Mesh.hpp"
#include "../params/ibvFunction.hpp"
#include "field.hpp"

/**
 * \brief Physical flux (allows to compute it at each node of the mesh).
 * \param fx Physical flux variable along x.
 * \param fy Physical flux variable along y.
 * \param C Parameter C of the Lax-Friedrichs numerical flux.
 * \param u Nodal vector of unknowns.
 */
void flux(Field& field);


/**
 * \brief Physical flux (allows to compute it at the boundary conditions).
 * \param fx Physical flux variable along x.
 * \param fy Physical flux variable along y.
 * \param u Nodal vector of unknowns.
 */
void flux(std::vector<double>& Fx, std::vector<double>& Fy, 
			const std::vector<double>& Q);


/**
 * \brief Function that allows to build the rhs of the DG method.
 * \param mesh The mesh of the problem.
 * \param I The rhs vector of the problem.
 * \param u Nodal vector of unknowns.
 * \param fx Physical flux variable along x, evaluated at each node.
 * \param fy Physical flux variable along y, evaluated at each node.
 * \param C Parameter C of the Lax-Friedrichs numerical flux.
 * \param factor Parameter that allows to get the weak or stong form.
 * \param t Time.
 * \param boundaries Information about the boundaries.
 */
void buildFlux(const Mesh& mesh, Eigen::VectorXd& IH,Eigen::VectorXd& IuH,
			Eigen::VectorXd& IvH, Field& field, double factor, double t, 
			const std::map<std::string, ibc>& boundaries);

#endif /* buildFlux_hpp */
