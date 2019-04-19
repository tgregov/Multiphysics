#include "phiPsi.hpp"


// see .hpp file for description
void LFShallow(const Edge& edge, Field& field, unsigned int j, double factor, 
                bool boundary, unsigned int indexJ, unsigned int indexFrontJ, 
                const SolverParams& solverParams)
{

    // gravity parameter
    double g = solverParams.fluxCoeffs[0];

    // computation of the absolute value of the eigenvalue lambda inside the element 
    double lambdaIn =   (field.u[1][indexJ]*edge.normal[0]
                    +   field.u[2][indexJ]*edge.normal[1])/field.u[0][indexJ];

    lambdaIn = (lambdaIn >= 0) ? lambdaIn + sqrt(g*field.u[0][indexJ]) :
                -lambdaIn + sqrt(g*field.u[0][indexJ]);

    // computation of the absolute value of the eigenvalue lambda outside the element 
    double lambdaOut;
    if(boundary)
    {
        lambdaOut = (field.uAtBC[1]*edge.normal[0] 
                        + field.uAtBC[2]*edge.normal[1])/field.uAtBC[0];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*field.uAtBC[0]) :
                    -lambdaOut + sqrt(g*field.uAtBC[0]);
    }
    else
    {
        lambdaOut = (field.u[1][indexFrontJ]*edge.normal[0]
                        + field.u[2][indexFrontJ]*edge.normal[1])
                    /field.u[0][indexFrontJ];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*field.u[0][indexFrontJ]) :
            -lambdaOut + sqrt(g*field.u[0][indexFrontJ]);
    }

    // computation of the value of C in the LF scheme
    double C = (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);

    // compute the numerical flux
    if(boundary)
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ] 
                        + field.FluxAtBC[dim][unk])/2
                    - C*edge.normal[dim]*(field.u[unk][indexJ] 
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
                    -(factor*field.flux[dim][unk][indexJ] 
                        + field.flux[dim][unk][indexFrontJ])/2
                    - C*edge.normal[dim]*(field.u[unk][indexJ] 
                        - field.u[unk][indexFrontJ])/2;
            }
        }
    }
}


// see .hpp file for description
void Roe(const Edge& edge, Field& field, unsigned int j, double factor, 
            bool boundary, unsigned int indexJ, unsigned int indexFrontJ, 
            const SolverParams& solverParams)
{

    // gravity parameter
    double g = solverParams.fluxCoeffs[0];

    // compute the Roe averages
    double hJSqrt = sqrt(field.u[0][indexJ]);
    double hFrontJSqrt = sqrt(field.u[0][indexFrontJ]);

    double uRoe = ((field.u[1][indexJ]/field.u[0][indexJ])*hJSqrt 
                    + (field.u[1][indexFrontJ]/field.u[0][indexFrontJ])*hFrontJSqrt)
                    /(hJSqrt + hFrontJSqrt);

    double vRoe =   ((field.u[2][indexJ]/field.u[0][indexJ])*hJSqrt
                    + (field.u[2][indexFrontJ]/field.u[0][indexFrontJ])*hFrontJSqrt)
                    /(hJSqrt + hFrontJSqrt);

    double cRoe = sqrt(g*(field.u[0][indexJ] + field.u[0][indexFrontJ])/2);

    // compute the (limited) Froude number
    double Fr = (uRoe*edge.normal[0] + vRoe*edge.normal[1])/cRoe;

    if(Fr < -1.0)
    {
        Fr = -1.0;
    }
    else if(Fr > 1.0)
    {
        Fr = 1.0;
    }

    // compute the numerical flux
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


// see .hpp file for description
void LFTransport(const Edge& edge, Field& field, unsigned int j, double factor, 
                    bool boundary, unsigned int indexJ, unsigned int indexFrontJ, 
                    const SolverParams& solverParams)
{

    // compute the value of the C of a pure transport LF scheme
    double C = fabs(solverParams.fluxCoeffs[0]*edge.normal[0] 
                    + solverParams.fluxCoeffs[1]*edge.normal[1]);

    // compute the numerical flux
    if(boundary)
    {
        for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
            {
                field.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ] 
                        + field.FluxAtBC[dim][unk])/2
                    - C*edge.normal[dim]*(field.u[unk][indexJ] 
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
                    -(factor*field.flux[dim][unk][indexJ] 
                        + field.flux[dim][unk][indexFrontJ])/2
                    - C*edge.normal[dim]*(field.u[unk][indexJ] 
                        - field.u[unk][indexFrontJ])/2;
            }
        }
    }
}


// see .hpp file for description
void mean(const Edge& edge, Field& field, unsigned int j, double factor, 
            bool boundary, unsigned int indexJ, unsigned int indexFrontJ, 
            const SolverParams& solverParams)
{

    // compute the numerical flux
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
