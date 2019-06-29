#ifndef transport_phiPsi_hpp_included
#define transport_phiPsi_hpp_included

#include "../../mesh/Mesh.hpp"
#include "../../params/Params.hpp"
#include "../../solver/field.hpp"


/**
 * \brief Function that computes the numerical LF flux for a pure transport case.
 * \param edge Edge of the current node (useful for the normal & the offset index).
 * \param field Structure containing the current unknowns of the DG-FEM.
 * \param partialField Structure containing temporary unknowns
 * (like here, the fluxes at the boundary).
 * \param j Index of the current node, with respect to the current element.
 * \param factor Parameter that determines the weak (+1) or strong form (-1).
 * \param boundary Boolean that specifies if we consider a boundary (1) or not (0).
 * \param indexJ Index of the current node, with respect to the whole mesh.
 * \param indexFrontJ Index of the oppsoite node, with respect to the whole mesh.
 * \param solverParams Structure containing the solver's parameters.
 */
void LFTransport(const Edge& edge, Field& field, PartialField& partialField,
					unsigned int j, double factor, bool boundary,
					unsigned int indexJ, unsigned int indexFrontJ,
                    const SolverParams& solverParams);

#endif // transport_phiPsi_hpp_included
