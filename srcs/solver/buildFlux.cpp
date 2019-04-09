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
	if(lambdaIn >= 0)
	{
		lambdaIn += sqrt(9.81*field.H[indexJ]);
	}
	else
	{
		lambdaIn = -lambdaIn + sqrt(9.81*field.H[indexJ]);
	}

	double lambdaOut = 	field.uH[indexFrontJ]/field.H[indexFrontJ]*normal[0]
					+	field.vH[indexFrontJ]/field.H[indexFrontJ]*normal[1];
	if(lambdaOut >= 0)
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
void buildFlux(const Mesh& mesh, Eigen::VectorXd& IH,Eigen::VectorXd& IuH,
			Eigen::VectorXd& IvH, Field& field, double factor, double t, 
			const std::map<std::string, ibc>& boundaries)
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
			Eigen::VectorXd partialIH(elmPropHD.nSF); partialIH.setZero();
			Eigen::VectorXd partialIuH(elmPropHD.nSF); partialIuH.setZero();
			Eigen::VectorXd partialIvH(elmPropHD.nSF); partialIvH.setZero(); // necessary

			// loop over the edges for the current element
			unsigned int nSigma = element.edges.size();
			for(unsigned int s = 0 ; s < nSigma ; ++s)
			{
				// current edge
				Edge edge = element.edges[s];

				// we first compute the matrix-vector product of dM with gx and gy
				Eigen::VectorXd gxH(elmPropHD.nSF); gxH.setZero();
				Eigen::VectorXd gyH(elmPropHD.nSF); gyH.setZero();
				Eigen::VectorXd dMgxH(elmPropHD.nSF), dMgyH(elmPropHD.nSF);

				Eigen::VectorXd gxuH(elmPropHD.nSF); gxuH.setZero();
				Eigen::VectorXd gyuH(elmPropHD.nSF); gyuH.setZero();
				Eigen::VectorXd dMgxuH(elmPropHD.nSF), dMgyuH(elmPropHD.nSF);

				Eigen::VectorXd gxvH(elmPropHD.nSF); gxvH.setZero();
				Eigen::VectorXd gyvH(elmPropHD.nSF); gyvH.setZero();
				Eigen::VectorXd dMgxvH(elmPropHD.nSF), dMgyvH(elmPropHD.nSF);

				for(unsigned int j = 0 ; j < edge.offsetInElm.size() ; ++j)
				{

					// global index of the current node
					unsigned int indexJ = element.offsetInU + edge.offsetInElm[j];

					// case of a boundary condition
					if (edge.edgeInFront.first == -1)
					{/*

						double fxAtBC, fyAtBC, uAtBC;
						ibc boundary = boundaries.at(edge.bcName);

						// node "in front" (at boundary condition)
						uAtBC = boundary.ibcFunc(edge.nodeCoordinate[j],
                               u[indexJ], t, boundary.coefficients);

						// physical flux "in front" (at boundary condition)
                        flux(fxAtBC, fyAtBC, uAtBC, fluxCoeffs);*/

                        // compute the numerical flux
                        // the weak/strong form is stored in "factor"
						gxH[edge.offsetInElm[j]] += 0;
						gxuH[edge.offsetInElm[j]] += 0;
						gxvH[edge.offsetInElm[j]] += 0;
						//-(factor*fx[indexJ] + fxAtBC)/2
						//	- C*edge.normal[0]*(u[indexJ] - uAtBC)/2;
						gyH[edge.offsetInElm[j]] += 0;
						gyuH[edge.offsetInElm[j]] += 0;
						gyvH[edge.offsetInElm[j]] += 0;
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
						gxH[edge.offsetInElm[j]] +=
							-(factor*field.FxH[indexJ] + field.FxH[indexFrontJ])/2
							- C*edge.normal[0]*(field.H[indexJ] - field.H[indexFrontJ])/2;
						gyH[edge.offsetInElm[j]] +=
							-(factor*field.FyH[indexJ] + field.FyH[indexFrontJ])/2
							- C*edge.normal[1]*(field.H[indexJ] - field.H[indexFrontJ])/2;

						gxuH[edge.offsetInElm[j]] +=
							-(factor*field.FxuH[indexJ] + field.FxuH[indexFrontJ])/2
							- C*edge.normal[0]*(field.uH[indexJ] - field.uH[indexFrontJ])/2;
						gyuH[edge.offsetInElm[j]] +=
							-(factor*field.FyuH[indexJ] + field.FyuH[indexFrontJ])/2
							- C*edge.normal[1]*(field.uH[indexJ] - field.uH[indexFrontJ])/2;

						gxvH[edge.offsetInElm[j]] +=
							-(factor*field.FxvH[indexJ] + field.FxvH[indexFrontJ])/2
							- C*edge.normal[0]*(field.vH[indexJ] - field.vH[indexFrontJ])/2;
						gyvH[edge.offsetInElm[j]] +=
							-(factor*field.FyvH[indexJ] + field.FyvH[indexFrontJ])/2
							- C*edge.normal[1]*(field.vH[indexJ] - field.vH[indexFrontJ])/2;
					}
				}

				// matrix-vector products between dM and gx/dy
				dMgxH = element.dM[s]*gxH;
				dMgyH = element.dM[s]*gyH;

				dMgxuH = element.dM[s]*gxuH;
				dMgyuH = element.dM[s]*gyuH;

				dMgxvH = element.dM[s]*gxvH;
				dMgyvH = element.dM[s]*gyvH;

				// dot product between dM and the normal
				partialIH += edge.determinantLD[0]*(
					edge.normal[0]*dMgxH + edge.normal[1]*dMgyH);

				partialIuH += edge.determinantLD[0]*(
					edge.normal[0]*dMgxuH + edge.normal[1]*dMgyuH);

				partialIvH += edge.determinantLD[0]*(
					edge.normal[0]*dMgxvH + edge.normal[1]*dMgyvH);				

			}

			// add the local rhs vector to the global one
			// [TO DO]: find some Eigen function that allows to do that efficiently
			for(unsigned int j = 0 ; j < elmPropHD.nSF ; ++j)
			{
				IH[element.offsetInU + j] = partialIH[j];
				IuH[element.offsetInU + j] = partialIuH[j];
				IvH[element.offsetInU + j] = partialIvH[j];
			}
		}
	}
}
