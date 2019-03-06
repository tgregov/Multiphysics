#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
#include <Eigen/Sparse>
#include "buildM.hpp"

/* Function that builds the [M] matrix of the DG method
 * Inputs:
 *  - meshParams: structure that contains all the geometric & mesh information
 *  - M: sparse matrix [M], defined before this function call
 *
 * Description:
 * The objective is to compute sum_k{w_k*l_i(x_k)*l_j(x_k)*det[J](x_k)}, where 
 * the sum is done over the Gauss points (GP). This is done in 3 steps:
 *  I. Precomputation of the elements T^(k)_{ij} = w_k*l_i(x_k)*l_j(x_k):
 *      Those elements being independent of the element, there are done once
 *      for all at the beginning. T^(k)_{ij} is stored as a vector of
 *      vector T, where the T[k] = T^(k) corresponds to one specific GP, 
 *      and is a vector containing the {ij} components. Since T^(k)_{ij} = 
 *      T^(k)_{ji} (i.e. this is symmetric), we only compute the upper-half
 *      part, such that if i and j goes from 0 to N-1, then there T[k] is
 *      a vector of length N*(N+1)/2, stored by starting with i = 0 and 
 *      going from j = 0 to N-1, then setting i = 1 and going from j = 1 to
 *      N-1, and so on. Furthemore, the {i,j} indices are stored for each 
 *      component T^(k)_{ij}, since it will be useful at the end to reconstruct
 *      the overal matrix.
 *  II. Computation of [M]_{ij} for each element:
 *      By looping over the elements, the matrices E^(e)_{ij}, which corresponds
 *      to the part of the overall matrix [M]_{ij} linked to the element e (i.e, 
 *      E^(e)_{ij} corresponds a block that appears in the diagonal of 
 *      [M]_{ij}), are determined. Given an element e, E^(e)_{ij} is calculated
 *      using a loop over the {ij} indices. Since those components are stored
 *      a 1D vector, we simply loop over l, where l goes along the vector
 *      components. Finally, the sum_k part of the initial formula is computed.
 *  III. Storage of the result using Eigen:
 *      The matrices E^(e)_{ij} are assembled and stored using a sparse format.
 */     

void buildM(const MeshParams& meshParams, Eigen::SparseMatrix<double>& M)
{

    // I. Precomputation of the elements T^(k)_{ij} = w_k*l_i(x_k)*l_j(x_k)
    std::vector<std::vector<double>> T;
    std::vector<std::pair<unsigned int, unsigned int>> IJ;

    for(unsigned int k = 0; k < meshParams.nGP; k++)
    {
        std::vector<double> t;

        for(unsigned int i = 0; i < meshParams.nSF; i++)
        {
            for(unsigned int j = i; j < meshParams.nSF; j++)
            {
                if(k == 0)
                {
                    IJ.push_back(std::pair<unsigned int, unsigned int>(i, j));
                }

                t.push_back(meshParams.intPoints[4*k + 3]
                    *meshParams.basisFunc[meshParams.nSF*k + i]
                    *meshParams.basisFunc[meshParams.nSF*k + j]);
            }
        }

        T.push_back(t);
    }

    // II. Computation of [M]_{ij} for each element
    std::vector<std::vector<double>> E;

    for(unsigned int elm = 0; elm < meshParams.nE; elm++)
    {
        std::vector<double> e;

        for(unsigned int l = 0; l < meshParams.nSF*(meshParams.nSF+1)/2; l++)
        {
            double sum = 0;
            for(unsigned int k = 0; k < meshParams.nGP; k++)
            {
                sum += T[k][l]*meshParams.determinant[elm*meshParams.nGP + k];
            }

            e.push_back(sum);
        }

        E.push_back(e);
    }

    // III. Storage of the result using Eigen
    std::vector<Eigen::Triplet<double>> index;

    for(unsigned int elm = 0; elm < meshParams.nE; elm++)
    {
        for(unsigned int l = 0; l < meshParams.nSF*(meshParams.nSF+1)/2; l++)
        {
            index.push_back(Eigen::Triplet<double>
                (IJ[l].first + elm*meshParams.nSF, 
                    IJ[l].second + elm*meshParams.nSF, 
                    E[elm][l]));

            if(IJ[l].first != IJ[l].second)
            {
                index.push_back(Eigen::Triplet<double>
                    (IJ[l].second + elm*meshParams.nSF, 
                        IJ[l].first + elm*meshParams.nSF, 
                        E[elm][l]));
            }
        }
    }

    M.setFromTriplets(index.begin(), index.end());
}
