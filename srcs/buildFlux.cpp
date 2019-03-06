#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI

#include <Eigen/Dense>
#include "readMesh.hpp"
#include <cmath>

void flux(Eigen::VectorXd<double>& fx, Eigen::VectorXd<double>& fy, double& C,
			const Eigen::VectorXd<double>& u)
{
	// first basic flux: simple transport
	std::vector<double> a;
	a.push_back(1.0);
	a.push_back(0.0);

	fx = a[0]*u;
	fy = a[1]*u;

	C = sqrt(a[0]*a[0] + a[1]*a[1]);
}


bool buildFlux(const MeshParams& meshParams, Eigen::VectorXd<double>& I,
				const Eigen::VectorXd<double>& u, const std::string& typeForm)
{


	Eigen::VectorXd<double> fx(meshParams.nE*meshParams.nSF), fy(meshParams.nE*meshParams.nSF);
	double C;
	flux(fx, fy, C, u);

	// loop over the elements
	for(unsigned int elm = 0; elm < meshParams.nE; elm++)
	{

		// we will build the vector I_j, j = 0, ..., nSF-1, corresponding to the
		// fluxes in the current element
		Eigen::VectorXd partialI(meshParams.nSF);
		partialI.setZero();

		std::vector<unsigned int> index = meshParams.index[elm];


		// loop, for each element, over the edges
		for(unsigned int s = 0; s < meshParams.nSigma; s++)
		{

			// we first compute the matrix-vector product of dM with gx and gy
			Eigen::VectorXd dMgx(meshParams.nSF), dMgy(meshParams.nSF);
			Eigen::VectorXd gx(meshParams.nSF), gy(meshParams.nSF);

			for(unsigned int j = 0; j < meshParams.nSF; j++)
			{

				// the type of form implies a different rhs vector
				float factor;

				if(typeForm.compare("strong"))
				{
					factor = -1.0;
				} else if(typeForm.compare("weak")){
					factor = +1.0;
				}
				else{
					std::cerr 	<< "The form  " << typeForm  << "does not exist !"
								<< std::endl;
					return false;
				}
				
				//DO NOT forget BC !!! 
				gx[j] = -(factor*fx[index[j]] + fx[opp(s, index[j])])/2
						- C*meshParams.normals[elm][s][0]*(u[index[j]] - u[meshParams.indexInFront(index[j])])/2;

				gy[j] = -(factor*fy[index[j]] + fy[opp(s, index[j])])/2
						- C*meshParams.normals[elm][s][1]*(u[index[j]] - u[meshParams.indexInFront(index[j]])/2;

			}

			dMgx = meshParams.dM[s]*gx;
			dMgy = meshParams.dM[s]*gy;

			// then we apply a scalar product and sum the current contribution
			// "+=" seems to work
			// /!\ nDet not disponible yet /!\
			partialI += meshParams.nDet[elm][s].first*dMgx + meshParams.nDet[elm][s].second*dMgy;

		}

		//Building of the vector I from the partialI
		// maybe there is some Eigen function that allows to do that
		for(unsigned int j = 0; j < meshParams.nSF; j++)
		{
			I[elm*meshParams.nSF + j] = partialI[j];
		}
	}

	return true;
}
