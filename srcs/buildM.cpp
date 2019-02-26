#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI
#include <gmsh.h>
#include <Eigen/Sparse>

bool buildM(Eigen::SparseMatrix<double>& M, const std::string& fileName, 
			const std::string& intScheme, const std::string& basisFunc)
{


 	

    // We need to compute:
    // sum_k{w_k*l_i(x_k)*l_j(x_k)*[detJ](x_k)}

    // parameters
    unsigned int nGP = intpts.size()/4; //maybe short
    unsigned int nSF = bf.size()/nGP;
    unsigned int nE = det.size()/nGP;

    // std::cout << "number of GP: " << nGP << std::endl; 
    // std::cout << "number of SF: " << nSF << std::endl; 
    // std::cout << "number of E: " << nE << std::endl; 

    // T^(k)_{ij}: k is the GP, {ij} is the l_i*l_j
    std::vector<std::vector<double>> T;

    for(unsigned int k = 0; k < nGP; k++)
    {
        std::vector<double> t;

        for(unsigned int i = 0; i < nSF; i++)
        {
            for(unsigned int j = i; j < nSF; j++)
            {
               t.push_back(intpts[4*k + 3]*bf[nGP*i + k]*bf[nGP*j + k]);
            }
        }

        T.push_back(t);
    }

    // E^(e)_{ij}: e is the element, {ij} is the l_i*l_j
    std::vector<std::vector<double>> E;

    for(unsigned int elm = 0; elm < nE; elm++)
    {
        std::vector<double> e;

        for(unsigned int l = 0; l < nSF*(nSF+1)/2; l++)
        {
            double sum = 0;
            for(unsigned int k = 0; k < nGP; k++)
            {
                sum += T[k][l]*det[elm*nGP + k];
            }

            e.push_back(sum);
        }

        E.push_back(e);
    }

    // assembly of the matrix M_ij
    // Eigen::SparseMatrix<double> M(nE*nSF, nE*nSF);
    std::vector<Eigen::Triplet<double>> index;

    for(unsigned int elm = 0; elm < nE; elm++)
    {   
        unsigned int i = 0;
        unsigned int j = 0;
        for(unsigned int l = 0; l < nSF*(nSF+1)/2; l++)
        {

            index.push_back(Eigen::Triplet<double>(i + elm*nSF, j + elm*nSF,
                                E[elm][l]));

            if(i != j)
            {
                index.push_back(Eigen::Triplet<double>(j + elm*nSF, i + elm*nSF,
                                E[elm][l]));
            }

            j += 1;
            if(j == nSF)
            {
                i += 1;
                j = i;
            }
        }
    }  

    M.setFromTriplets(index.begin(), index.end());
    std::cout << M << std::endl;

	return true;
}