#include <iostream>
#include <cmath>
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


bool buildFlux(const Mesh2D& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u,
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy, const double& C,
	const std::string& typeForm, unsigned int numNodes)
{

	// the type of form is stored in factor
	float factor;
	if(!typeForm.compare("strong")) 	factor = -1.0;
	else if(!typeForm.compare("weak"))	factor = +1.0;
	else
	{
		std::cerr 	<< "The form  " << typeForm  << "does not exist !"
					<< std::endl;
		return false;
	}

	// loop over the entities
	for(unsigned int ent = 0 ; ent < mesh.entities.size() ; ent++)
	{
		// current entity
		Entity2D entity = mesh.entities[ent];

		// loop over the elements
		for(unsigned int elm = 0 ; elm < entity.elements.size() ; elm++)
		{
			// current element
			Element2D element = entity.elements[elm];

			// get the properties of the current element type
            ElementProperty elmProp1D = mesh.elementProperties1D.at(element.elementType1D);
            ElementProperty elmProp2D = mesh.elementProperties2D.at(element.elementType2D);

			// partial rhs vector
			Eigen::VectorXd partialI(elmProp2D.nSF); partialI.setZero();


			// I. BUILD THE DELTAM MATRIX
			// compute sum_k{w_k*l_i*l_j}
			// [TO DO]: only works for linear SF
			// [TO DO]: put it in mesh2D
			double lala = 0.0;
			double lalb = 0.0;
			for (unsigned int k = 0 ; k < elmProp1D.nGP ; ++k)
			{
				lala += elmProp1D.prodFunc[k][0];
				lalb += elmProp1D.prodFunc[k][1];
			}

			// compute the indices of the components
			unsigned int nSigma = element.edges.size();
			std::vector<Eigen::SparseMatrix<double>> dM;


			for(unsigned int s = 0 ; s < nSigma ; s++)
			{
				Eigen::SparseMatrix<double> dMs(elmProp2D.nSF, elmProp2D.nSF);
				std::vector<Eigen::Triplet<double>> indices;

				indices.push_back(Eigen::Triplet<double>(s, s, lala));
				indices.push_back(Eigen::Triplet<double>(s, (s+1) % nSigma, lalb));
				indices.push_back(Eigen::Triplet<double>((s+1) % nSigma, s, lalb));
				indices.push_back(Eigen::Triplet<double>((s+1) % nSigma, (s+1) % nSigma, lala));
				dMs.setFromTriplets(indices.begin(), indices.end());

				dM.push_back(dMs);
			}


			// II. COMPUTE THE RHS
			// [TO DO] optmize x2
			// loop, for each element, over the edges
			for(unsigned int s = 0 ; s < nSigma ; ++s)
			{
				// current edge
				Edge edge = element.edges[s];

				// we first compute the matrix-vector product of dM with gx and gy
				Eigen::VectorXd gx(elmProp2D.nSF), gy(elmProp2D.nSF);
				Eigen::VectorXd dMgx(elmProp2D.nSF), dMgy(elmProp2D.nSF);

				// loop over the SF
				for(unsigned int j = 0 ; j < elmProp2D.nSF ; ++j)
				{
                    // [TO DO]: neighbours in other entities
                    unsigned int frontOffsetInU = entity.elements[std::get<0>(edge.edgeInFront)].offsetInU;
                    unsigned int frontJ = std::get<1>(edge.edgeInFront);

					// [TO DO]: only works for linear elements
					// Boundary condition case
					if (std::get<0>(edge.edgeInFront) == -1)
					{
						if(elm == 1)
						{
                            gx[j] = 1;

                            gy[j] = 0;
						}
						else
						{
                            gx[j] = -(factor*fx[element.offsetInU + j] + 0)/2
                                    - C*element.edges[s].normal.first*(u[element.offsetInU + j] - 0)/2;

                            gy[j] = 0;
						}
					}
					else
					{
						gx[j] = -(factor*fx[element.offsetInU + j] + fx[frontOffsetInU + frontJ])/2
								- C*element.edges[s].normal.first*(u[element.offsetInU + j] - u[frontOffsetInU + frontJ])/2;

						gy[j] = -(factor*fy[element.offsetInU + j] + fy[frontOffsetInU + frontJ])/2
								- C*element.edges[s].normal.second*(u[element.offsetInU + j] - u[frontOffsetInU + frontJ])/2;
					}
				}

				dMgx = dM[s]*gx;
				dMgy = dM[s]*gy;

				// then we apply a scalar product and sum the current contribution
				// "+=" seems to work
				// [TO DO]: constant determinant
				partialI += edge.determinant1D[0]*(
					element.edges[s].normal.first*dMgx
					+ element.edges[s].normal.second*dMgy);
			}

			// Building of the vector I from the partialI
			// maybe there is some Eigen function that allows to do that
			for(unsigned int j = 0 ; j < elmProp2D.nSF ; ++j)
			{
				I[element.offsetInU + j] = partialI[j];
				std::cout << "elm = " << elm << " | node = " << j << " => " << partialI[j] << std::endl;
			}
		}
	}

	return true;
}
