#include <iostream>
#include "flux.hpp"

void fluxShallow(Field& field, const SolverParams& solverParams, bool boundary)
{
    double g = solverParams.fluxCoeffs[0];

    if(boundary)
    {
        field.FluxAtBC[0][0] = field.uAtBC[1];
        field.FluxAtBC[1][0] = field.uAtBC[2];

        field.FluxAtBC[0][1] = field.uAtBC[1]*field.uAtBC[1]/field.uAtBC[0] + g/2*field.uAtBC[0]*field.uAtBC[0];
        field.FluxAtBC[0][2] = field.uAtBC[1]*field.uAtBC[2]/field.uAtBC[0];

        field.FluxAtBC[1][1] = field.uAtBC[1]*field.uAtBC[2]/field.uAtBC[0];
        field.FluxAtBC[1][2] = field.uAtBC[2]*field.uAtBC[2]/field.uAtBC[0] + g/2*field.uAtBC[0]*field.uAtBC[0];
    }
    else
    {
        if (field.u[0].minCoeff() <= 0)
        {
            std::cerr << "WARNING: Negative height or division by 0 !" << std::endl;
        }

        field.flux[0][0] = field.u[1];
        field.flux[1][0] = field.u[2];

        field.flux[0][1]  = field.u[1].array().square()/field.u[0].array()
                    + g/2*field.u[0].array().square();
        field.flux[0][2]  = field.u[1].array()*field.u[2].array()/field.u[0].array();

        field.flux[1][1]  = field.u[1].array()*field.u[2].array()/field.u[0].array();
        field.flux[1][2]  = field.u[2].array().square()/field.u[0].array()
                    + g/2*field.u[0].array().square();
    }
}

