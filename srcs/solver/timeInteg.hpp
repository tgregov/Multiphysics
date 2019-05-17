#ifndef timeInteg_hpp_included
#define timeInteg_hpp_included

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../mesh/Mesh.hpp"
#include "../params/Params.hpp"


/**
 * \brief Time integrate the equations (DG-FEM)
 * \param mesh The mesh representing the domain of interest
 * \param solverParams The structure in which the parameters of the solver are.
 * \param fileName The name of the .msh file containing the mesh.
 * \param resultsName name of the .msh file that will contain the results
 * \return true if time integration happened without problems, false otherwise.
 */
bool timeInteg(const Mesh& mesh, SolverParams& solverParams,
				const std::string& fileName, const std::string& resultsName);

#endif /* timeInteg_hpp */
