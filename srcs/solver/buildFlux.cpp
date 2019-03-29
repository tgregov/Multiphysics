#include <iostream>
#include <cmath>
#include "buildFlux.hpp"


// see .hpp file for description
void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, const Eigen::VectorXd& u,
          const std::vector<double>& fluxCoeffs)
{
	// first basic flux: simple transport
	fx = fluxCoeffs[0]*u;
	fy = fluxCoeffs[1]*u;
}


// see .hpp file for description
void flux(double& fx, double& fy, double u, const std::vector<double>& fluxCoeffs)
{
	// first basic flux: simple transport
    fx = fluxCoeffs[0]*u;
	fy = fluxCoeffs[1]*u;
}

static double computeC(const std::vector<double>& normal,
                       const std::vector<double>& fluxCoeffs)
{
    return fabs(fluxCoeffs[0]*normal[0] + fluxCoeffs[1]*normal[1]);
}


// see .hpp file for description
void buildFlux(const Mesh& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u,
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy,
	double factor, double t, const std::map<std::string, ibc>& boundaries,
	const std::vector<double>& fluxCoeffs)
{

	// loop over the entities
	for(size_t ent = 0 ; ent < mesh.entities.size() ; ent++)
	{
		// current entity
		Entity entity = mesh.entities[ent];

		// loop over the elements
		for(size_t elm = 0 ; elm < entity.elements.size() ; elm++)
		{
			// current element
			Element element = entity.elements[elm];

			// get the properties of the current element type
            ElementProperty elmPropLD
            	= mesh.elementProperties.at(element.elementTypeLD);
            ElementProperty elmPropHD
            	= mesh.elementProperties.at(element.elementTypeHD);

			// local I vector for the current element
			Eigen::VectorXd partialI(elmPropHD.nSF); partialI.setZero(); // necessary

			// loop over the edges for the current element
			unsigned int nSigma = element.edges.size();
			for(unsigned int s = 0 ; s < nSigma ; ++s)
			{
				// current edge
				Edge edge = element.edges[s];

				// we first compute the matrix-vector product of dM with gx and gy
				Eigen::VectorXd gx(elmPropHD.nSF); gx.setZero();
				Eigen::VectorXd gy(elmPropHD.nSF); gy.setZero();
				Eigen::VectorXd dMgx(elmPropHD.nSF), dMgy(elmPropHD.nSF);


				double C = computeC(edge.normal, fluxCoeffs);

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
						uAtBC = boundary.ibcFunc(edge.nodeCoordinate[j][0],
							edge.nodeCoordinate[j][1],
							0.0, u[indexJ], t, boundary.coefficients);

						// physical flux "in front" (at boundary condition)
                        flux(fxAtBC, fyAtBC, uAtBC, fluxCoeffs);

                        // compute the numerical flux
                        // the weak/strong form is stored in "factor"
						gx[edge.offsetInElm[j]] += -(factor*fx[indexJ] + fxAtBC)/2
							- C*edge.normal[0]*(u[indexJ] - uAtBC)/2;
						gy[edge.offsetInElm[j]] += -(factor*fy[indexJ] + fyAtBC)/2
							- C*edge.normal[1]*(u[indexJ] - uAtBC)/2;
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
							- C*edge.normal[0]*(u[indexJ] - u[indexFrontJ])/2;
						gy[edge.offsetInElm[j]] +=
							-(factor*fy[indexJ] + fy[indexFrontJ])/2
							- C*edge.normal[1]*(u[indexJ] - u[indexFrontJ])/2;
					}
				}

				// matrix-vector products between dM and gx/dy
				dMgx = element.dM[s]*gx;
				dMgy = element.dM[s]*gy;

				// dot product between dM and the normal
				partialI += edge.determinantLD[0]*(
					edge.normal[0]*dMgx + edge.normal[1]*dMgy);
			}

			// add the local rhs vector to the global one
			// [TO DO]: find some Eigen function that allows to do that efficiently
			for(unsigned int j = 0 ; j < elmPropHD.nSF ; ++j)
			{
				I[element.offsetInU + j] = partialI[j];
			}
		}
	}
}
