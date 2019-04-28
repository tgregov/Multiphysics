#include <iostream>
#include "source.hpp"

void sourceShallowCst(Field& field, const SolverParams& solverParams)
{
    field.s[0].setZero();
    field.s[1] = solverParams.sourceCoeffs[0]*field.u[2]
    + solverParams.fluxCoeffs[0]*solverParams.sourceCoeffs[1]*field.u[0];
    field.s[2] = -solverParams.sourceCoeffs[0]*field.u[1];
    + solverParams.fluxCoeffs[0]*solverParams.sourceCoeffs[2]*field.u[0];
}

void sourceShallowLinCst(Field& field, const SolverParams& solverParams)
{
    field.s[0].setZero();
    field.s[1] = solverParams.sourceCoeffs[0]*solverParams.fluxCoeffs[1]*field.u[2];
    field.s[2] = -solverParams.sourceCoeffs[0]*solverParams.fluxCoeffs[1]*field.u[1];
}
