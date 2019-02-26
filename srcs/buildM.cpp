#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
#include <gmsh.h>
#include "buildM.hpp"

void buildM(const MeshParams& meshParams, Eigen::SparseMatrix<double>& M)
{
    // T^(k)_{ij}:
    // k is the GP
    // {ij} is the l_i*l_j
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

                t.push_back(meshParams.intPoints[4*k + 3]*meshParams.basisFunc[meshParams.nGP*i + k]*meshParams.basisFunc[meshParams.nGP*j + k]);
            }
        }

        T.push_back(t);
    }

    // E^(e)_{ij}
    // e is the element
    // {ij} is the l_i*l_j
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

    // assembly of the matrix M_ij
    std::vector<Eigen::Triplet<double>> index;

    for(unsigned int elm = 0; elm < meshParams.nE; elm++)
    {
        for(unsigned int l = 0; l < meshParams.nSF*(meshParams.nSF+1)/2; l++)
        {
            index.push_back(Eigen::Triplet<double>(IJ[l].first + elm*meshParams.nSF, IJ[l].second + elm*meshParams.nSF, E[elm][l]));

            if(IJ[l].first != IJ[l].second)
            {
                index.push_back(Eigen::Triplet<double>(IJ[l].second + elm*meshParams.nSF, IJ[l].first + elm*meshParams.nSF, E[elm][l]));
            }
        }
    }

    M.setFromTriplets(index.begin(), index.end());
}
