#ifndef timeInteg_hpp_included
#define timeInteg_hpp_included

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "mesh/Mesh2D.hpp"
#include "Solver.hpp"

/**
 * \brief Time integrate the equations (DG-FEM)
 * \param mesh The mesh representing the domain of interest
 * \param solverParams The structure in which the parameters of the solver are.
 * \param fileName The name of the file containing the mesh.
 * \return true if time integration happened without problems, false otherwise.
 */
bool timeInteg(const Mesh2D& mesh, const SolverParams& solverParams,
	             const std::string& fileName);

#endif /* timeInteg_hpp */
