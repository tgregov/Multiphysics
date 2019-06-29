#ifndef transport_flux_hpp_included
#define transport_flux_hpp_included

#include "../../params/Params.hpp"
#include "../../solver/field.hpp"


/**
 * \brief Function that computes the physical flux for a pure transport.
 * \param field Structure containing all the information about the computed unknowns.
 * \param partialField Structure containing temporary unknowns.
 * \param solverParams Structure containing the solver's parameters.
 * \param boundary Boolean that specifies if we consider a boundary (1) or not (0).
 */
void fluxTransport(Field& field, PartialField& partialField,
					const SolverParams& solverParams, bool boundary);

#endif // transport_flux_hpp_included
