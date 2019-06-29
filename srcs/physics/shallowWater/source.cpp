#include "source.hpp"


// see .hpp file for description
void sourceShallowCstGradCstFrict(Field& field, const SolverParams& solverParams)
{
    field.s[0].setZero();

    field.s[1] = solverParams.sourceCoeffs[1]*field.u[2]
            + solverParams.fluxCoeffs[0]*solverParams.sourceCoeffs[2]*field.u[0]
            - solverParams.sourceCoeffs[4]*field.u[1]/solverParams.sourceCoeffs[0];

    field.s[2] = -solverParams.sourceCoeffs[1]*field.u[1]
            + solverParams.fluxCoeffs[0]*solverParams.sourceCoeffs[3]*field.u[0]
            - solverParams.sourceCoeffs[4]*field.u[2]/solverParams.sourceCoeffs[0];
}


// see .hpp file for description
void sourceShallowCstGradQuadFrict(Field& field, const SolverParams& solverParams)
{
    Eigen::VectorXd speedNorm(field.u[0].size());
    speedNorm = field.u[1].array()*field.u[1].array()
                + field.u[2].array()*field.u[2].array();
    speedNorm = speedNorm.array().sqrt();
    speedNorm = speedNorm.array()/field.u[0].array();

    field.s[0].setZero();

    field.s[1] = solverParams.sourceCoeffs[1]*field.u[2].array()
        + solverParams.fluxCoeffs[0]*solverParams.sourceCoeffs[2]*field.u[0].array()
        - solverParams.sourceCoeffs[4]*speedNorm.array()*field.u[1].array()
        /(solverParams.sourceCoeffs[0]*field.u[0].array());

    field.s[2] = -solverParams.sourceCoeffs[1]*field.u[1].array()
        + solverParams.fluxCoeffs[0]*solverParams.sourceCoeffs[3]*field.u[0].array()
        - solverParams.sourceCoeffs[4]*speedNorm.array()*field.u[2].array()
        /(solverParams.sourceCoeffs[0]*field.u[0].array());
}
