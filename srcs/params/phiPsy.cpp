#include "phiPsy.hpp"

void LF(const Edge& edge, Field& field, unsigned int j,
        double factor, bool boundary, unsigned int indexJ,
        unsigned int indexFrontJ, const SolverParams& solverParams)
{

    // compute the C constant of LF
    double g = solverParams.fluxCoeffs[0];

    double lambdaIn =   (field.u[1][indexJ]*edge.normal[0]
                    +   field.u[2][indexJ]*edge.normal[1])/field.u[0][indexJ];

    lambdaIn = (lambdaIn >= 0) ? lambdaIn + sqrt(g*field.u[0][indexJ]) :
                -lambdaIn + sqrt(g*field.u[0][indexJ]);

    double lambdaOut;
    if(boundary)
    {
        lambdaOut =     (field.uAtBC[1]*edge.normal[0] + field.uAtBC[2]*edge.normal[1])/field.uAtBC[0];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*field.uAtBC[0]) :
                    -lambdaOut + sqrt(g*field.uAtBC[0]);
    }
    else
    {
        lambdaOut =     (field.u[1][indexFrontJ]*edge.normal[0]
                + field.u[2][indexFrontJ]*edge.normal[1])/field.u[0][indexFrontJ];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*field.u[0][indexFrontJ]) :
            -lambdaOut + sqrt(g*field.u[0][indexFrontJ]);
    }

    double C = (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);

    // compute the flux
    if(boundary)
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ] + field.FluxAtBC[dim][unk])/2
                    - C*edge.normal[dim]*(field.u[unk][indexJ] - field.uAtBC[unk])/2;
            }
        }
    }
    else
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ] + field.flux[dim][unk][indexFrontJ])/2
                    - C*edge.normal[dim]*(field.u[unk][indexJ] - field.u[unk][indexFrontJ])/2;
            }
        }
    }
}


void Roe(const Edge& edge, Field& field, unsigned int j,
         double factor, bool boundary, unsigned int indexJ,
         unsigned int indexFrontJ, const SolverParams& solverParams)
{
    // compute the Roe averages
    double hJSqrt = sqrt(field.u[0][indexJ]);
    double hFrontJSqrt = sqrt(field.u[0][indexFrontJ]);

    double uRoe =   ((field.u[1][indexJ]/field.u[0][indexJ])*hJSqrt
                    +
                    (field.u[1][indexFrontJ]/field.u[0][indexFrontJ])*hFrontJSqrt)
                    /(hJSqrt + hFrontJSqrt);

    double vRoe =   ((field.u[2][indexJ]/field.u[0][indexJ])*hJSqrt
                    +
                    (field.u[2][indexFrontJ]/field.u[0][indexFrontJ])*hFrontJSqrt)
                    /(hJSqrt + hFrontJSqrt);

    double g = solverParams.fluxCoeffs[0];
    double cRoe = sqrt(g*(field.u[0][indexJ] + field.u[0][indexFrontJ])/2);

    double Fr = (uRoe*edge.normal[0] + vRoe*edge.normal[1])/cRoe;
    // std::cout << "\nFr = " << Fr << std::endl;

    if(Fr < -1.0)
    {
        Fr = -1.0;
    }
    else if(Fr > 1.0)
    {
        Fr = 1.0;
    }

    // compute the flux
    if(boundary)
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                     -((Fr + factor)*field.flux[dim][unk][indexJ]+
                     (1 - Fr)*field.FluxAtBC[dim][unk])/2
                     - cRoe*(1-Fr*Fr)*edge.normal[dim]*(field.u[unk][indexJ]
                     - field.uAtBC[unk])/2;
            }
        }
    }
    else
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                     -((Fr + factor)*field.flux[dim][unk][indexJ]+
                     (1 - Fr)*field.flux[dim][unk][indexFrontJ])/2
                     - cRoe*(1-Fr*Fr)*edge.normal[dim]*(field.u[unk][indexJ]
                     - field.u[unk][indexFrontJ])/2;
            }
        }
    }
}


void mean(const Edge& edge, Field& field, unsigned int j,
          double factor, bool boundary, unsigned int indexJ,
          unsigned int indexFrontJ, const SolverParams& solverParams)
{

    // compute the flux
    if(boundary)
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ]
                    + field.FluxAtBC[dim][unk])/2;
            }
        }
    }
    else
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ]
                    + field.flux[dim][unk][indexFrontJ])/2;
            }
        }
    }
}
