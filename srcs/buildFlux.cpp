#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Mesh2D.hpp"

void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, double& C,
			const Eigen::VectorXd& u)
{
	// first basic flux: simple transport
	std::vector<double> a;
	a.push_back(1.0);
	a.push_back(0.0);

	fx = a[0]*u;
	fy = a[1]*u;

	C = sqrt(a[0]*a[0] + a[1]*a[1]);
}


bool buildFlux(Mesh2D& mesh, Eigen::VectorXd& I,
				const Eigen::VectorXd& u, const std::string& typeForm, unsigned int numNodes)
{
	//loop over the entities
	for(unsigned int ent = 0 ; ent < mesh.entities.size() ; ent++)
	{
		Entity2D entity = mesh.entities[ent];
		Eigen::VectorXd fx(numNodes), fy(numNodes);
		double C;
		flux(fx, fy, C, u);

		//loop over the elements
		for(unsigned int elm = 0 ; elm < entity.elements.size() ; elm++)
		{
			Element2D element = entity.elements[elm];

			// get the properties of the current element type
            ElementProperty elmProp = mesh.elementProperties1D[element.elementType1D];
            ElementProperty elmProp2D = mesh.elementProperties2D[element.elementType2D];


			// we will build the vector I_j, j = 0, ..., nSF-1, corresponding to the
			// fluxes in the current element
			Eigen::VectorXd partialI(elmProp2D.nSF);
			partialI.setZero();



			
			double lala = 0.0;
			double lalb = 0.0;
			//only works with linear shape functions
			for (unsigned int k = 0 ; k < elmProp.nGP ; ++k)
			{
				lala += elmProp.pondFunc[k][0]*elmProp.pondFunc[k][0];
				lalb += elmProp.pondFunc[k][1]*elmProp.pondFunc[k][0];
			}
		

			//Computation of matrix delta M

			//creation of matrix M
			Eigen::SparseMatrix<double> dM(elmProp2D.nSF, elmProp2D.nSF);
			std::vector<Eigen::Triplet<double>> indices;

			unsigned int nSigma = element.edges.size();
			for(unsigned int s = 0; s < nSigma; s++)
			{
				indices.push_back(Eigen::Triplet<double>(s, s, lala));
				indices.push_back(Eigen::Triplet<double>(s, (s+1) % nSigma, lalb));
				indices.push_back(Eigen::Triplet<double>((s+1) % nSigma, s, lalb));
				indices.push_back(Eigen::Triplet<double>((s+1) % nSigma, (s+1) % nSigma, lala));
			}
			dM.setFromTriplets(indices.begin(), indices.end());


			
			// loop, for each element, over the edges
			
			for(unsigned int s = 0 ; s < nSigma ; ++s)
			{


				
				Edge edge = element.edges[s];
				
				// we first compute the matrix-vector product of dM with gx and gy
				Eigen::VectorXd dMgx(elmProp2D.nSF), dMgy(elmProp2D.nSF);
				Eigen::VectorXd gx(elmProp2D.nSF), gy(elmProp2D.nSF);

				

				
				for(unsigned int j = 0 ; j < elmProp2D.nSF ; ++j)
				{

					// the type of form implies a different rhs vector
					float factor;

					if(typeForm.compare("strong"))
					{
						factor = -1.0;
					} else if(typeForm.compare("weak"))
					{
						factor = +1.0;
					}
					else{
						std::cerr 	<< "The form  " << typeForm  << "does not exist !"
									<< std::endl;
						return false;
					}
					// /!\Check for nonlinear elements
					if (std::get<0>(edge.edgeInFront) == -1)
					{
						gx[j] = 0.0;
						gy[j] = 0.0;
					}else{
						unsigned int frontOffsetInU = entity.elements[std::get<0>(edge.edgeInFront)].offsetInU;
						unsigned int frontJ = std::get<1>(edge.edgeInFront);
						//DO NOT forget BC !!!

						gx[j] = -(factor*fx(u[element.offsetInU + j]) + fx(u[frontOffsetInU + frontJ]))/2
								- C*element.edgesNormal[s].first*(u[element.offsetInU + j]) - u[frontOffsetInU + frontJ]/2;

						gy[j] = -(factor*fy(u[element.offsetInU + j]) + fy(u[frontOffsetInU + frontJ]))/2
								- C*element.edgesNormal[s].second*u[element.offsetInU + j] - u[frontOffsetInU + frontJ]/2;
					}

				}
				
				
				dMgx = dM*gx;
				dMgy = dM*gy;

				// then we apply a scalar product and sum the current contribution
				// "+=" seems to work
				// constant determinant 
				partialI += edge.determinant1D[0]*(element.edgesNormal[s].first*dMgx + element.edgesNormal[s].second*dMgx);

			}

			//Building of the vector I from the partialI
			// maybe there is some Eigen function that allows to do that
			
			for(unsigned int j = 0 ; j < elmProp2D.nSF ; ++j)
			{
				I[element.offsetInU + j] = partialI[j];
			}

		}
	}

	return true;
}
