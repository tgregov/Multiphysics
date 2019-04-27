#include <iostream>
#include "source.hpp"

void sourceShallow(Field& field, const SolverParams& solverParams)
{
    field.s[0].setZero();
    field.s[1] = solverParams.sourceCoeffs[0]*field.u[2];
    field.s[2] = -solverParams.sourceCoeffs[0]*field.u[1];
}

void sourceShallowLin(Field& field, const SolverParams& solverParams)
{
    field.s[0].setZero();
    field.s[1] = solverParams.sourceCoeffs[0]*solverParams.fluxCoeffs[1]*field.u[2];
    field.s[2] = -solverParams.sourceCoeffs[0]*solverParams.fluxCoeffs[1]*field.u[1];
}
