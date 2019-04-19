#ifndef flux_hpp_included
#define flux_hpp_included

#include "Params.hpp"
#include "../solver/field.hpp"


void fluxShallow(Field& field, const SolverParams& solverParams, bool boundary);
void fluxTransport(Field& field, const SolverParams& solverParams, bool boundary);

#endif // flux_hpp_included
