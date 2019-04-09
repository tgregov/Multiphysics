#include <iostream>
#include <cmath>
#include "buildFlux.hpp"

// see .hpp file for description
void flux(Field& field)
{

	// TO DO: add Joachim's solution
	// TO DO: div by zero warning
	for(size_t i = 0 ; i < field.H.size() ; ++i)
	{
		field.FxH[i] = field.uH[i];
		field.FyH[i] = field.vH[i];

		field.FxuH[i] = (field.uH[i]*field.uH[i])/field.H[i] 
						+ 9.81/2*field.H[i]*field.H[i];
		field.FxvH[i] = field.uH[i]*field.vH[i]/field.H[i];

		field.FyuH[i] = field.uH[i]*field.vH[i]/field.H[i];
		field.FyvH[i] = (field.vH[i]*field.vH[i])/field.H[i] 
						+ 9.81/2*field.H[i]*field.H[i];
	}
}


// see .hpp file for description
void flux(std::vector<double>& Fx, std::vector<double>& Fy, 
			const std::vector<double>& Q)
{
	Fx[0] = Q[1];
	Fy[0] = Q[2];
	
	Fx[1] = Q[1]*Q[1]/Q[0] + 9.81/2*Q[0]*Q[0];
	Fx[2] = Q[1]*Q[2]/Q[0]; 

	Fy[1] = Q[1]*Q[2]/Q[0];
	Fy[2] = Q[2]*Q[2]/Q[0] + 9.81/2*Q[0]*Q[0];
	 
}

static double computeC(const std::vector<double>& normal,
                       const Field& field, 
                       unsigned int indexJ,
                       unsigned int indexFrontJ)
{

	double lambdaIn = 	field.uH[indexJ]/field.H[indexJ]*normal[0]
					+	field.vH[indexJ]/field.H[indexJ]*normal[1];
	if(lambdaIn > = 0)
	{
		lambdaIn += sqrt(9.81*field.H[indexJ]);
	}
	else
	{
		lambdaIn = -lambdaIn + sqrt(9.81*field.H[indexJ]);
	}

	double lambdaOut = 	field.uH[indexFrontJ]/field.H[indexFrontJ]*normal[0]
					+	field.vH[indexFrontJ]/field.H[indexFrontJ]*normal[1];
	if(lambdaOut > = 0)
	{
		lambdaOut += sqrt(9.81*field.H[indexFrontJ]);
	}
	else
	{
		lambdaOut = -lambdaOut + sqrt(9.81*field.H[indexFrontJ]);
	}

    return (lambdaIn > lambdaOut ? lambdaIn : lambdaOut);
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
						uAtBC = boundary.ibcFunc(edge.nodeCoordinate[j],
                               u[indexJ], t, boundary.coefficients);

						// physical flux "in front" (at boundary condition)
                        flux(fxAtBC, fyAtBC, uAtBC, fluxCoeffs);

                        // compute the numerical flux
                        // the weak/strong form is stored in "factor"
						gx[edge.offsetInElm[j]] += 0;
						//-(factor*fx[indexJ] + fxAtBC)/2
						//	- C*edge.normal[0]*(u[indexJ] - uAtBC)/2;
						gy[edge.offsetInElm[j]] += 0;
						//-(factor*fy[indexJ] + fyAtBC)/2
						//	- C*edge.normal[1]*(u[indexJ] - uAtBC)/2;
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

						double C = computeC(edge.normal, field, indexJ, indexFrontJ);

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
