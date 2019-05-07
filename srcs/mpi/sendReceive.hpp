#ifndef sendReceive_hpp_included
#define sendReceive_hpp_included

#include "../solver/field.hpp"
#include "../utils/utils.hpp"
#include "../params/Params.hpp"
#include "../mesh/Mesh.hpp"


void exchangeFlux(const Field& field, CompleteField& compField,
                  const DomainDiv& domainDiv, unsigned int rank,
                  const SolverParams& solverParams, const Mesh& mesh);

void exchangeUnk(const Field& field, CompleteField& compField,
                  const DomainDiv& domainDiv, unsigned int rank,
                  const SolverParams& solverParams, const Mesh& mesh);

#endif /* sendReceive_hpp_included */
