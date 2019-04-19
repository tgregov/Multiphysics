#ifndef flux_hpp_included
#define flux_hpp_included

#include "../params/Params.hpp"
#include "../solver/field.hpp"


/**
 * \brief Function that computes the physical flux for the shallow water equations.
 * \param field Structure containing all the information about the computed unknowns.
 * \param solverParams Structure containing the solver's parameters.
 * \param boundary Boolean that specifies if we consider a boundary (1) or not (0).
 */
void fluxShallow(Field& field, const SolverParams& solverParams, bool boundary);


/**
 * \brief Function that computes the physical flux for a pure transport.
 * \param field Structure containing all the information about the computed unknowns.
 * \param solverParams Structure containing the solver's parameters.
 * \param boundary Boolean that specifies if we consider a boundary (1) or not (0).
 */
void fluxTransport(Field& field, const SolverParams& solverParams, bool boundary);

#endif // flux_hpp_included
