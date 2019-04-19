#include <iostream>
#include <cmath>
#include "buildFlux.hpp"
#include <Eigen/Dense>

// see .hpp file for description
void buildFlux(const Mesh& mesh, Field& field, double factor, double t,
               const SolverParams& solverParams)
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
			for(unsigned short unk = 0 ; unk < field.partialIu.size() ; ++unk)
            {
                field.partialIu[unk].resize(elmPropHD.nSF);
                field.partialIu[unk].setZero();
            }

			// loop over the edges for the current element
			unsigned int nSigma = element.edges.size();
			for(unsigned int s = 0 ; s < nSigma ; ++s)
			{
				// current edge
				Edge edge = element.edges[s];

				// we first compute the matrix-vector product of dM with gx and gy
				for(unsigned short dim = 0 ; dim < field.g.size() ; ++dim)
                {
                    for(unsigned short unk = 0 ; unk < field.g[dim].size() ; ++unk)
                    {
                            field.g[dim][unk].resize(elmPropHD.nSF);
                            field.g[dim][unk].setZero();
                    }
                }

				for(unsigned int j = 0 ; j < edge.offsetInElm.size() ; ++j)
				{
					// global index of the current node
					unsigned int indexJ = element.offsetInU + edge.offsetInElm[j];

					// case of a boundary condition
					if (edge.edgeInFront.first == -1)
					{
					    //[TO DO]: Fix boundary for shallow water
					    for(unsigned short unk = 0 ; unk < field.uForBC.size() ; ++unk)
                            field.uForBC[unk] = field.u[unk][indexJ];

						ibc boundary = solverParams.boundaryConditions.at(edge.bcName);
						boundary.ibcFunc(field.uAtBC, edge.nodeCoordinate[j], t,
                                         field.uForBC, edge.normal, boundary.coefficients);

                        solverParams.flux(field, solverParams, true);

                        // compute the numerical flux
                        // the weak/strong form is stored in "factor"
                        solverParams.phiPsy(edge, field, j, factor,
                                            true, indexJ, 0, solverParams);
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
                        solverParams.phiPsy(edge, field, j, factor,
                                            false, indexJ, indexFrontJ, solverParams);
					}
				}

				// dot product between dM and the normal
				for(unsigned short unk = 0 ; unk < field.partialIu.size() ; ++unk)
                {
                    field.partialIu[unk]+= edge.determinantLD[0]*(
					edge.normal[0]*element.dM[s]*field.g[0][unk] + edge.normal[1]*element.dM[s]*field.g[1][unk]);
                }
			}

			// add the local rhs vector to the global one
			// [TO DO]: find some Eigen function that allows to do that efficiently
			for(unsigned short unk = 0 ; unk < field.Iu.size() ; ++unk)
            {
                for(unsigned int j = 0 ; j < elmPropHD.nSF ; ++j)
                {
                    field.Iu[unk][element.offsetInU + j] = field.partialIu[unk][j];
                }
            }
		}
	}
}
