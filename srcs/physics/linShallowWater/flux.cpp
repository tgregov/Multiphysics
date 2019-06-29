#include <iostream>
#include "flux.hpp"

// see .hpp file for description
void fluxShallowLin(Field& field, PartialField& partialField,
                    const SolverParams& solverParams, bool boundary)
{

    // the physical flux is given by
    //  Fx = [h0*u,     g*h0*H, 0       ]
    //  Fy = [h0*v,     0,      g*h_0*h ]

    if(boundary)
    {
        // flux for the (x, y) coordinates for H
        partialField.FluxAtBC[0][0] = partialField.uAtBC[1];
        partialField.FluxAtBC[1][0] = partialField.uAtBC[2];

        // flux for the (x, y) coordinates for H*u
        partialField.FluxAtBC[0][1] = solverParams.fluxCoeffs[0]
                                *partialField.uAtBC[0]*solverParams.fluxCoeffs[1];
        partialField.FluxAtBC[0][2] = 0.0;

        // flux for the (x, y) coordinates for H*v
        partialField.FluxAtBC[1][1] = 0.0;
        partialField.FluxAtBC[1][2] = solverParams.fluxCoeffs[0]
                                *partialField.uAtBC[0]*solverParams.fluxCoeffs[1];
    }
    else
    {
        if (field.u[0].minCoeff() <= 0)
        {
            std::cerr << "WARNING: Negative height or division by 0 !" << std::endl;
        }

        // flux for the (x, y) coordinates for H
        field.flux[0][0] = field.u[1];
        field.flux[1][0] = field.u[2];

        // flux for the (x, y) coordinates for u
        field.flux[0][1] = solverParams.fluxCoeffs[0]*solverParams.fluxCoeffs[1]
                                *field.u[0];
        // field.flux[0][2] = Eigen::VectorXd::Zero(field.u[0].size());

        // flux for the (x, y) coordinates for v
        // field.flux[1][1] = Eigen::VectorXd::Zero(field.u[0].size());
        field.flux[1][2] = solverParams.fluxCoeffs[0]*solverParams.fluxCoeffs[1]
                                *field.u[0];
    }
}
