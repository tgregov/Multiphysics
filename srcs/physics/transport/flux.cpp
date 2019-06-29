#include <iostream>
#include "flux.hpp"


// see .hpp file for description
void fluxTransport(Field& field, PartialField& partialField,
                    const SolverParams& solverParams, bool boundary)
{
    // the physical flux is F = a*Q", where "a" is a 2D vector
    if(boundary)
    {
        // flux for the x coordinate
        partialField.FluxAtBC[0][0] = solverParams.fluxCoeffs[0]
                                        *partialField.uAtBC[0];
        // flux for the y coordinate
        partialField.FluxAtBC[1][0] = solverParams.fluxCoeffs[1]
                                        *partialField.uAtBC[0];
    }
    else
    {
        // flux for the x coordinate
        field.flux[0][0] = solverParams.fluxCoeffs[0]*field.u[0];
        // flux for the y coordinate
        field.flux[1][0] = solverParams.fluxCoeffs[1]*field.u[0];
    }
}
