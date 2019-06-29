#include "meanPhiPsi.hpp"

void mean(const Edge& edge, Field& field, PartialField& partialField, unsigned int j, double factor,
            bool boundary, unsigned int indexJ, unsigned int indexFrontJ,
            const SolverParams& solverParams)
{

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
                        + partialField.FluxAtBC[dim][unk])/2;
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
                        + field.flux[dim][unk][indexFrontJ])/2;
            }
        }
    }
}
