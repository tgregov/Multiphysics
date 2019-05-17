#ifndef buildFlux_hpp
#define buildFlux_hpp

#include "../mesh/Mesh.hpp"
#include "../params/Params.hpp"
#include "../solver/field.hpp"


/**
 * \brief Function that allows to build the rhs of the DG method.
 * \param mesh The mesh of the problem.
 * \param field Structure containing all the information about the computed unknowns.
 * \param factor Parameter for a strong or weak form of DG method (+1 or -1).
 * \param t Current time of the simulation.
 * \param solverParams Structure containing the solver's parameters.
 */
void buildFlux(const Mesh& mesh, Field& field, double factor, double t,
               const SolverParams& solverParams);

#endif /* buildFlux_hpp */
