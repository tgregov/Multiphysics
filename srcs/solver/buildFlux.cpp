#include <iostream>
#include <cmath>
#include "buildFlux.hpp"


// see .hpp file for description
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


// see .hpp file for description
void flux(double& fx, double& fy, double u)
{
	// first basic flux: simple transport
	double ax = 1.0;
	double ay = 0.0;

	fx = ax*u;
	fy = ay*u;
}


// see .hpp file for description
void buildFlux(const Mesh2D& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u,
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy, double C,
	double factor, unsigned int numNodes, double t,
	const std::map<std::string, ibc>& boundaries)
{

	// loop over the entities
	for(size_t ent = 0 ; ent < mesh.entities.size() ; ent++)
	{
		// current entity
		Entity2D entity = mesh.entities[ent];

		// loop over the elements
		for(size_t elm = 0 ; elm < entity.elements.size() ; elm++)
		{
			// current element
			Element2D element = entity.elements[elm];

			// get the properties of the current element type
            ElementProperty elmProp1D
            	= mesh.elementProperties1D.at(element.elementType1D);
            ElementProperty elmProp2D
            	= mesh.elementProperties2D.at(element.elementType2D);

			// local I vector for the current element
			Eigen::VectorXd partialI(elmProp2D.nSF); partialI.setZero(); // necessary

			// loo over the edges for the current element
			unsigned int nSigma = element.edges.size();
			for(unsigned int s = 0 ; s < nSigma ; ++s)
			{
				// current edge
				Edge edge = element.edges[s];

				// we first compute the matrix-vector product of dM with gx and gy
				Eigen::VectorXd gx(elmProp2D.nSF); gx.setZero();
				Eigen::VectorXd gy(elmProp2D.nSF); gy.setZero();
				Eigen::VectorXd dMgx(elmProp2D.nSF), dMgy(elmProp2D.nSF);

				for(unsigned int j = 0 ; j < edge.offsetInElm.size() ; ++j)
				{

					// global index of the current node
					unsigned int indexJ = element.offsetInU + edge.offsetInElm[j];

					// case of a boundary condition
					if (edge.edgeInFront.first == -1)
					{

						double fxAtBC, fyAtBC, uAtBC;
						ibc boundary = boundaries.at(edge.bcName);

						// node "in front" (at boundary condition)
						uAtBC = boundary.ibcFunc(edge.nodeCoordinate[j].first,
							edge.nodeCoordinate[j].second,
							0.0, u[indexJ], t, boundary.coefficients);

						// physical flux "in front" (at boundary condition)
                        flux(fxAtBC, fyAtBC, uAtBC);

                        // compute the numerical flux
                        // the weak/strong form is stored in "factor"
						gx[edge.offsetInElm[j]] += -(factor*fx[indexJ] + fxAtBC)/2
							- C*edge.normal.first*(u[indexJ] - uAtBC)/2;
						gy[edge.offsetInElm[j]] += -(factor*fy[indexJ] + fyAtBC)/2
							- C*edge.normal.second*(u[indexJ] - uAtBC)/2;
					}
					else // general case
					{
					    // [TO DO]: neighbours in other entities
					    // global index of the node "in front"
                       	unsigned int indexFrontJ =
                       				entity
                       					.elements[edge.edgeInFront.first]
                       					.offsetInU
                       				+ entity
                       					.elements[edge.edgeInFront.first]
                       					.edges[edge.edgeInFront.second]
                       					.offsetInElm[edge.nodeIndexEdgeInFront[j]];

                        // compute the numerical flux
                        // the weak/strong form is stored in "factor"
						gx[edge.offsetInElm[j]] +=
							-(factor*fx[indexJ] + fx[indexFrontJ])/2
							- C*edge.normal.first*(u[indexJ] - u[indexFrontJ])/2;
						gy[edge.offsetInElm[j]] +=
							-(factor*fy[indexJ] + fy[indexFrontJ])/2
							- C*edge.normal.second*(u[indexJ] - u[indexFrontJ])/2;
					}
				}

				// matrix-vector products between dM and gx/dy
				dMgx = element.dM[s]*gx;
				dMgy = element.dM[s]*gy;

				// dot product between dM and the normal
				partialI += edge.determinant1D[0]*(
					edge.normal.first*dMgx + edge.normal.second*dMgy);
			}

			// add the local rhs vector to the global one
			// [TO DO]: find some Eigen function that allows to do that efficiently
			for(unsigned int j = 0 ; j < elmProp2D.nSF ; ++j)
			{
				I[element.offsetInU + j] = partialI[j];
			}
		}
	}
}
