#include "phiPsy.hpp"

double LF(const std::vector<double>& edgeNormal, const Field& field,
          unsigned int dim, unsigned int unk, double factor, bool boundary,
          unsigned int indexJ, unsigned int indexFrontJ, double C)
{
    if(boundary)
    {
        return -(factor*field.flux[dim][unk][indexJ] + field.FluxAtBC[dim][unk])/2
               - C*edgeNormal[dim]*(field.u[unk][indexJ] - field.uAtBC[unk])/2;
    }
    else
    {
        return -(factor*field.flux[dim][unk][indexJ] + field.flux[dim][unk][indexFrontJ])/2
               - C*edgeNormal[dim]*(field.u[unk][indexJ] - field.u[unk][indexFrontJ])/2;
    }
}

double LFShallowC(const std::vector<double>& edgeNormal, const Field& field, bool boundary,
                  unsigned int indexJ, unsigned int indexFrontJ, const SolverParams& solverParams)
{
    double g = solverParams.fluxCoeffs[0];

    double lambdaIn = 	(field.u[1][indexJ]*edgeNormal[0]
					+	field.u[2][indexJ]*edgeNormal[1])/field.u[0][indexJ];

    lambdaIn = (lambdaIn >= 0) ? lambdaIn + sqrt(g*field.u[0][indexJ]) :
                -lambdaIn + sqrt(g*field.u[0][indexJ]);

    double lambdaOut;
    if(boundary)
    {
        lambdaOut = 	(field.uAtBC[1]*edgeNormal[0] + field.uAtBC[2]*edgeNormal[1])/field.uAtBC[0];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*field.uAtBC[0]) :
                    -lambdaOut + sqrt(g*field.uAtBC[0]);
    }
    else
    {
        lambdaOut = 	(field.u[1][indexFrontJ]*edgeNormal[0]
                + field.u[2][indexFrontJ]*edgeNormal[1])/field.u[0][indexFrontJ];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*field.u[0][indexFrontJ]) :
            -lambdaOut + sqrt(g*field.u[0][indexFrontJ]);
    }

    return (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);
}

double meanC(const std::vector<double>& edgeNormal, const Field& field, bool boundary,
             unsigned int indexJ, unsigned int indexFrontJ, const SolverParams& solverParams)
{
    return 0.0;
}
