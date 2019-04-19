#ifndef phiPsy_hpp_included
#define phiPsy_hpp_included

#include "../mesh/Mesh.hpp"
#include "Params.hpp"
#include "../solver/field.hpp"

void LF(const Edge& edge, Field& field, unsigned int j,
        double factor, bool boundary, unsigned int indexJ,
        unsigned int indexFrontJ, const SolverParams& solverParams);

void Roe(const Edge& edgel, Field& field, unsigned int j,
         double factor, bool boundary, unsigned int indexJ,
         unsigned int indexFrontJ, const SolverParams& solverParams);

void mean(const Edge& edge, Field& field, unsigned int j,
            double factor, bool boundary, unsigned int indexJ,
            unsigned int indexFrontJ, const SolverParams& solverParams);
#endif // phiPsy_hpp_included
