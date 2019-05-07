#include "phiPsi.hpp"


// see .hpp file for description
void LFShallow(const Edge& edge, Field& field, PartialField& partialField,
               const CompleteField& compField, unsigned int j, double factor,
               bool boundary, unsigned int indexJ, unsigned int indexFrontJ,
               const SolverParams& solverParams, unsigned int nodePrec)
{
    // gravity parameter
    double g = solverParams.fluxCoeffs[0];

    // computation of the absolute value of the eigenvalue lambda inside the element
    double lambdaIn =   (field.u[1][indexJ - nodePrec]*edge.normal[0]
                    +   field.u[2][indexJ - nodePrec]*edge.normal[1])/field.u[0][indexJ - nodePrec];

    lambdaIn = (lambdaIn >= 0) ? lambdaIn + sqrt(g*field.u[0][indexJ - nodePrec]) :
                -lambdaIn + sqrt(g*field.u[0][indexJ - nodePrec]);

    double lambdaOut;

    // compute the numerical flux
    if(boundary)
    {
        // computation of the absolute value of the eigenvalue lambda outside the element
        lambdaOut = (partialField.uAtBC[1]*edge.normal[0]
                        + partialField.uAtBC[2]*edge.normal[1])/partialField.uAtBC[0];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*partialField.uAtBC[0]) :
                    -lambdaOut + sqrt(g*partialField.uAtBC[0]);

        // computation of the value of C in the LF scheme
        double C = (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);

        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ-nodePrec]
                        + partialField.FluxAtBC[dim][unk]
                    + C*edge.normal[dim]*(field.u[unk][indexJ-nodePrec]
                        - partialField.uAtBC[unk]))/2;
            }
        }
    }
    else
    {
        // computation of the absolute value of the eigenvalue lambda outside the element
        lambdaOut = (field.u[1][indexFrontJ]*edge.normal[0]
                        + field.u[2][indexFrontJ]*edge.normal[1])
                    /field.u[0][indexFrontJ];

        lambdaOut = (lambdaOut >= 0) ? lambdaOut + sqrt(g*field.u[0][indexFrontJ]) :
            -lambdaOut + sqrt(g*field.u[0][indexFrontJ]);

        // computation of the value of C in the LF scheme
        double C = (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);

        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ-nodePrec]
                        + compField.flux[dim][unk][indexFrontJ]
                    + C*edge.normal[dim]*(field.u[unk][indexJ-nodePrec]
                        - compField.u[unk][indexFrontJ]))/2;
            }
        }
    }
}


// see .hpp file for description
void Roe(const Edge& edge, Field& field, PartialField& partialField,
               const CompleteField& compField, unsigned int j, double factor,
               bool boundary, unsigned int indexJ, unsigned int indexFrontJ,
               const SolverParams& solverParams, unsigned int nodePrec)
{
    // gravity parameter
    double g = solverParams.fluxCoeffs[0];

    double hJSqrt, hFrontJSqrt, uRoe, vRoe, cRoe;

    // compute the numerical flux
    if(boundary)
    {
        // compute the Roe averages
        hJSqrt = sqrt(field.u[0][indexJ]);
        hFrontJSqrt = sqrt(partialField.uAtBC[0]);

        uRoe = ((field.u[1][indexJ]/field.u[0][indexJ])*hJSqrt
                        + (partialField.uAtBC[1]/partialField.uAtBC[0])*hFrontJSqrt)
                        /(hJSqrt + hFrontJSqrt);

        vRoe =   ((field.u[2][indexJ]/field.u[0][indexJ])*hJSqrt
                        + (partialField.uAtBC[2]/partialField.uAtBC[0])*hFrontJSqrt)
                        /(hJSqrt + hFrontJSqrt);

        cRoe = sqrt(g*(field.u[0][indexJ] + partialField.uAtBC[0])/2);

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

        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -((Fr + factor)*field.flux[dim][unk][indexJ]+
                        (1 - Fr)*partialField.FluxAtBC[dim][unk]
                    + cRoe*(1-Fr*Fr)*edge.normal[dim]*(field.u[unk][indexJ]
                        - partialField.uAtBC[unk]))/2;
            }
        }
    }
    else
    {
        // compute the Roe averages
        hJSqrt = sqrt(field.u[0][indexJ]);
        hFrontJSqrt = sqrt(field.u[0][indexFrontJ]);

        uRoe = ((field.u[1][indexJ]/field.u[0][indexJ])*hJSqrt
                        + (field.u[1][indexFrontJ]/field.u[0][indexFrontJ])*hFrontJSqrt)
                        /(hJSqrt + hFrontJSqrt);

        vRoe =   ((field.u[2][indexJ]/field.u[0][indexJ])*hJSqrt
                        + (field.u[2][indexFrontJ]/field.u[0][indexFrontJ])*hFrontJSqrt)
                        /(hJSqrt + hFrontJSqrt);

        cRoe = sqrt(g*(field.u[0][indexJ] + field.u[0][indexFrontJ])/2);

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

        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -((Fr + factor)*field.flux[dim][unk][indexJ]+
                        (1 - Fr)*field.flux[dim][unk][indexFrontJ]
                    + cRoe*(1-Fr*Fr)*edge.normal[dim]*(field.u[unk][indexJ]
                        - field.u[unk][indexFrontJ]))/2;
            }
        }
    }
}


// see .hpp file for description
void LFTransport(const Edge& edge, Field& field, PartialField& partialField,
               const CompleteField& compField, unsigned int j, double factor,
               bool boundary, unsigned int indexJ, unsigned int indexFrontJ,
               const SolverParams& solverParams, unsigned int nodePrec)
{

    // compute the value of the C of a pure transport LF scheme
    double C = fabs(solverParams.fluxCoeffs[0]*edge.normal[0]
                    + solverParams.fluxCoeffs[1]*edge.normal[1]);

    // compute the numerical flux
    if(boundary)
    {
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ - nodePrec]
                        + partialField.FluxAtBC[dim][unk]
                    + C*edge.normal[dim]*(field.u[unk][indexJ - nodePrec]
                        - partialField.uAtBC[unk]))/2;
            }
        }
    }
    else
    {
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ - nodePrec]
                        + field.flux[dim][unk][indexFrontJ]
                    + C*edge.normal[dim]*(field.u[unk][indexJ - nodePrec]
                        - field.u[unk][indexFrontJ]))/2;
            }
        }
    }
}


// see .hpp file for description
void mean(const Edge& edge, Field& field, PartialField& partialField,
               const CompleteField& compField, unsigned int j, double factor,
               bool boundary, unsigned int indexJ, unsigned int indexFrontJ,
               const SolverParams& solverParams, unsigned int nodePrec)
{

    // compute the numerical flux
    if(boundary)
    {
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ - nodePrec]
                        + partialField.FluxAtBC[dim][unk])/2;
            }
        }
    }
    else
    {
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ - nodePrec]
                        + field.flux[dim][unk][indexFrontJ])/2;
            }
        }
    }
}
