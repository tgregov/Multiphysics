#ifndef computeError_hpp_included
#define computeError_hpp_included

#include <string>
#include "../mesh/Mesh.hpp"
#include "../params/Params.hpp"


/**
 * \brief compute the L2 error of an approximate solution
 * \param mesh The mesh representing the domain of interest
 * \param solverParams the solver parameter
 * \param resultsName the .msh file that contains the results
 * \param meshName the .msh file containing the mesh
 */
bool computeError(const Mesh& mesh, const SolverParams& solverParams, const std::string& meshName,
				 const std::string& reslutsName, double& error);

#endif /* computeError_hpp */