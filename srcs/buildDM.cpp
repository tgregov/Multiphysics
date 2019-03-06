#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
#include <Eigen/Sparse>
#include "buildDM.hpp"


void buildDM(MeshParams& meshParams)
{
	double lala = 0.0;
	double lalb = 0.0;
	unsigned int nSigma = meshParams.nSigma;

	for (unsigned int k = 0 ; k < meshParams.nGPInf ; k++)
	{
	    lala += meshParams.intPointsInferior[4*k + 3]
	                    *meshParams.basisFuncInferior[meshParams.nSFInf*k]
	                    *meshParams.basisFuncInferior[meshParams.nSFInf*k];

		lalb += meshParams.intPointsInferior[4*k + 3]
	                    *meshParams.basisFuncInferior[meshParams.nSFInf*k]
	                    *meshParams.basisFuncInferior[meshParams.nSFInf*k + 1];
	}


	for(unsigned int s = 0 ; s < meshParams.nSigma ; s++)
	{
		meshParams.dM.push_back(Eigen::SparseMatrix<double>(meshParams.nSF,meshParams.nSF));

		std::vector<Eigen::Triplet<double>> index;
		index.push_back(Eigen::Triplet<double>(s, s, lala));
		index.push_back(Eigen::Triplet<double>(s, (s+1) % nSigma, lalb));
		index.push_back(Eigen::Triplet<double>((s+1) % nSigma, s, lalb));
		index.push_back(Eigen::Triplet<double>((s+1) % nSigma, (s+1)%nSigma, lala));

		meshParams.dM[s].setFromTriplets(index.begin(), index.end());

	}

}

