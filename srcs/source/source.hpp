#ifndef source_hpp_included
#define source_hpp_included

#include <Eigen/Dense>
#include "../params/Params.hpp"
#include "../solver/field.hpp"

void sourceShallowCst(Field& field, const SolverParams& solverParams);

void sourceShallowLinCst(Field& field, const SolverParams& solverParams);


#endif // source_hpp_included
