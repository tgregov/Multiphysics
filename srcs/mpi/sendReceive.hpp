#ifndef sendReceive_hpp_included
#define sendReceive_hpp_included

#include "../solver/field.hpp"
#include "../utils/utils.hpp"
#include "../params/Params.hpp"
#include "../mesh/Mesh.hpp"

/**
 * \brief Exchange the physical flux between the MPI threads.
 * \param field Structure containing all the information about the computed unknowns
 * (MPI thread local).
 * \param compField Structure containing the whole unknowns and fluxes vector.
 * \param domainDiv Structure representing how the nodes
 * are split into the MPI threads.
 * \param rank Rank of the MPI thread.
 * \param solverParams Structure containing the solver's parameters.
 * \param mesh The structure that contains the mesh.
 */
void exchangeFlux(const Field& field, CompleteField& compField,
                  const DomainDiv& domainDiv, unsigned int rank,
                  const SolverParams& solverParams, const Mesh& mesh);


/**
 * \brief Exchange the unknowns between the MPI threads.
 * \param field Structure containing all the information about the computed unknowns
 * (MPI thread local).
* \param compField Structure containing the whole unknowns and fluxes vector.
 * \param domainDiv Structure representing how the nodes
 * are split into the MPI threads.
 * \param rank Rank of the MPI thread.
 * \param solverParams Structure containing the solver's parameters.
 * \param mesh The structure that contains the mesh.
 */
void exchangeUnk(const Field& field, CompleteField& compField,
                  const DomainDiv& domainDiv, unsigned int rank,
                  const SolverParams& solverParams, const Mesh& mesh);

#endif /* sendReceive_hpp_included */
