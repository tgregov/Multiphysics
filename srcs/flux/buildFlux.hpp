#ifndef buildFlux_hpp
#define buildFlux_hpp

#include "../mesh/Mesh.hpp"
#include "../params/Params.hpp"
#include "../solver/field.hpp"
#include "../utils/utils.hpp"


/**
 * \brief Function that allows to build the rhs of the DG method.
 * \param mesh The mesh of the problem.
 * \param field Structure containing all the information about the computed unknowns
 * (MPI thread local).
 * \param compField Structure containing the whole unknowns and fluxes vector.
 * \param factor Parameter for a strong or weak form of DG method (+1 or -1).
 * \param t Current time of the simulation
 * \param solverParams Structure containing the solver's parameters.
 * \param domainDiv Structure representing how the nodes
 * are split into the MPI threads.
 * \param rank Rank of the MPI thread.
 */
void buildFlux(const Mesh& mesh, Field& field, const CompleteField& compField,
               double factor, double t, const SolverParams& solverParams,
               const DomainDiv& domainDiv, unsigned int rank);

#endif /* buildFlux_hpp */
