#include "phiPsi.hpp"


// see .hpp file for description
void LFTransport(const Edge& edge, Field& field, PartialField& partialField, unsigned int j, double factor,
                    bool boundary, unsigned int indexJ, unsigned int indexFrontJ,
                    const SolverParams& solverParams)
{

    // compute the value of the C of a pure transport LF scheme
    double C = fabs(solverParams.fluxCoeffs[0]*edge.normal[0]
                    + solverParams.fluxCoeffs[1]*edge.normal[1]);

    // compute the numerical flux
    if(boundary)
    {
        //computation of g
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ]
                        + partialField.FluxAtBC[dim][unk]
                    + C*edge.normal[dim]*(field.u[unk][indexJ]
                        - partialField.uAtBC[unk]))/2;
            }
        }
    }
    else
    {
        //computation of g
        for(unsigned short dim = 0 ; dim < partialField.g.size() ; ++dim)
        {
            for(unsigned short unk = 0 ; unk < partialField.g[dim].size() ; ++unk)
            {
                partialField.g[dim][unk][edge.offsetInElm[j]] +=
                    -(factor*field.flux[dim][unk][indexJ]
                        + field.flux[dim][unk][indexFrontJ]
                    + C*edge.normal[dim]*(field.u[unk][indexJ]
                        - field.u[unk][indexFrontJ]))/2;
            }
        }
    }
}
