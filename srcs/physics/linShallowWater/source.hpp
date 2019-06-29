#ifndef linShallow_source_hpp_included
#define linShallow_source_hpp_included

#include <Eigen/Dense>
#include "../../params/Params.hpp"
#include "../../solver/field.hpp"


/**
 * \brief Function that computes the source term of the shallow water equation
 * (coriolis).
 * \param field Structure containing all the information about the computed unknowns.
 * \param solverParams Structure containing the solver's parameters.
 */
void sourceShallowLinCst(Field& field, const SolverParams& solverParams);

#endif // linShallow_source_hpp_included
