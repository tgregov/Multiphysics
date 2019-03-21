#include <iostream>
#include <cmath>
#include <Eigen/Sparse>
#include "./mesh/Mesh2D.hpp"


void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, double& C,
			const Eigen::VectorXd& u)
{
	// first basic flux: simple transport
	double ax = 1.0;
	double ay = 0.0;

	fx = ax*u;
	fy = ay*u;

	C = sqrt(ax*ax + ay*ay);
}


void flux(double& fx, double& fy, double u)
{
	// first basic flux: simple transport
	double ax = 1.0;
	double ay = 0.0;

	fx = ax*u;
	fy = ay*u;
}


double valueAtBC(const std::string& bcName,
					double x, double y, double z, double t)
{
	if(!bcName.compare("BC_Left"))
	{
		return exp(-(t-2)*(t-2)/0.1);
	}
	else
	{
		return 0.0;
	}
}


bool buildFlux(const Mesh2D& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u,
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy, const double& C,
	double factor, unsigned int numNodes, double t)
{

	// the type of form is stored in factor

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
            ElementProperty elmProp1D 
            	= mesh.elementProperties1D.at(element.elementType1D);
            ElementProperty elmProp2D 
            	= mesh.elementProperties2D.at(element.elementType2D);

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
				indices.push_back(Eigen::Triplet<double>((s+1) % nSigma, 
									(s+1) % nSigma, lala));
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

				gx.setZero();
				gy.setZero();

				for(unsigned int j = 0 ; j < edge.offsetInElm.size() ; ++j)
				{

					unsigned int indexJ = element.offsetInU + edge.offsetInElm[j];
					if (edge.edgeInFront.first == -1)
					{
						double fxAtBC, fyAtBC, uAtBC;
						uAtBC = valueAtBC(edge.bcName,
							edge.nodeCoordinate[j].first,
							edge.nodeCoordinate[j].second,
							0.0, t);

                        flux(fxAtBC, fyAtBC, uAtBC);

						gx[edge.offsetInElm[j]] += -(factor*fx[indexJ] + fxAtBC)/2
							- C*element.edges[s].normal.first*(u[indexJ] - uAtBC)/2;
						gy[edge.offsetInElm[j]] += -(factor*fy[indexJ] + fyAtBC)/2
							- C*element.edges[s].normal.second*(u[indexJ] - uAtBC)/2;
					}
					else
					{
					    // [TO DO]: neighbours in other entities
                       	unsigned int indexFrontJ = entity
                       					.elements[edge.edgeInFront.first]
                       					.offsetInU
                       				+ entity
                       					.elements[edge.edgeInFront.first]
                       					.edges[edge.edgeInFront.second]
                       					.offsetInElm[edge.nodeIndexEdgeInFront[j]];

						gx[edge.offsetInElm[j]] += 
							-(factor*fx[indexJ] + fx[indexFrontJ])/2
							- C*element.edges[s].normal.first*(u[indexJ] 
																- u[indexFrontJ])/2;
						gy[edge.offsetInElm[j]] += 
							-(factor*fy[indexJ] + fy[indexFrontJ])/2
							- C*element.edges[s].normal.second*(u[indexJ] 
																- u[indexFrontJ])/2;
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
			}
		}
	}

	return true;
}
