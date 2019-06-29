#ifndef shallow_source_hpp_included
#define shallow_source_hpp_included

#include <Eigen/Dense>
#include "../../params/Params.hpp"
#include "../../solver/field.hpp"


/**
 * \brief Function that computes the source term of the shallow water equation
 * (linear gradient of bed and linear bottom friction term, plus coriolis).
 * \param field Structure containing all the information about the computed unknowns.
 * \param solverParams Structure containing the solver's parameters.
 */
 void sourceShallowCstGradCstFrict(Field& field, const SolverParams& solverParams);


/**
 * \brief Function that computes the source term of the shallow water equation
 * (linear gradient of bed and quadratic bottom friction term, plus coriolis).
 * \param field Structure containing all the information about the computed unknowns.
 * \param solverParams Structure containing the solver's parameters.
 */
void sourceShallowCstGradQuadFrict(Field& field, const SolverParams& solverParams);


#endif // shallow_source_hpp_included
