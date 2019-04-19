#ifndef phiPsy_hpp_included
#define phiPsy_hpp_included

#include <vector>
#include "Params.hpp"
#include "../solver/field.hpp"

double LF(const std::vector<double>& edgeNormal, const Field& field,
          unsigned int dim, unsigned int unk, double factor, bool boundary,
          unsigned int indexJ, unsigned int indexFrontJ, double C);

double LFShallowC(const std::vector<double>& edgeNormal, const Field& field, bool boundary,
                  unsigned int indexJ, unsigned int indexFrontJ, const SolverParams& solverParams);

double meanC(const std::vector<double>& edgeNormal, const Field& field, bool boundary,
             unsigned int indexJ, unsigned int indexFrontJ, const SolverParams& solverParams);

#endif // phiPsy_hpp_included
