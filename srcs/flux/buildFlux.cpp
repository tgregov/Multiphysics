#include <iostream>
#include <cmath>
#include "buildFlux.hpp"
#include <Eigen/Dense>

// see .hpp file for description
void buildFlux(const Mesh& mesh, Field& field, const CompleteField& compField,
               double factor, double t, const SolverParams& solverParams,
               const DomainDiv& domainDiv, unsigned int rank)
{
	// loop over the entities
	for(size_t ent = 0 ; ent < mesh.entities.size() ; ent++)
	{
	    Entity entity = mesh.entities[ent];
		// loop over the elements
		//#pragma omp parallel for default(none) shared(field, compField, mesh, entity, solverParams, factor, t, nodePrec)
		for(size_t elm = domainDiv.elementPrec[rank] ; elm < domainDiv.elementPrec[rank] + domainDiv.element[rank] ; elm++)
		{
			PartialField partialField(solverParams.nUnknowns, mesh.dim);

			// local I vector for the current element
			for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
            {
                partialField.partialIu[unk].resize(mesh.elementProperties.at(entity.elements[elm].elementTypeHD).nSF);
                partialField.partialIu[unk].setZero();
            }

			// loop over the edges for the current element
			unsigned int nSigma = entity.elements[elm].edges.size();
			for(unsigned int s = 0 ; s < nSigma ; ++s)
			{
				// current edge
				// we first compute the matrix-vector product of dM with gx and gy
				for(unsigned short dim = 0 ; dim < mesh.dim ; ++dim)
                {
                    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
                    {
                            partialField.g[dim][unk].resize(mesh.elementProperties.at(entity.elements[elm].elementTypeHD).nSF);
                            partialField.g[dim][unk].setZero();
                    }
                }

				for(unsigned int j = 0 ; j < entity.elements[elm].edges[s].offsetInElm.size() ; ++j)
				{
					// global index of the current node
					unsigned int indexJ = entity.elements[elm].offsetInU + entity.elements[elm].edges[s].offsetInElm[j];

					// case of a boundary condition
					if (entity.elements[elm].edges[s].edgeInFront.first == -1)
					{
                        // compute the boundary condition
						ibc boundary
							= solverParams.boundaryConditions.at(entity.elements[elm].edges[s].bcName);
						boundary.ibcFunc(partialField.uAtBC, entity.elements[elm].edges[s].nodeCoordinate[j], t,
											field, indexJ-domainDiv.nodePrec[rank], entity.elements[elm].edges[s].normal,
											boundary.coefficients, solverParams.fluxCoeffs);

                        solverParams.flux(field, partialField, solverParams, true);

                        // compute the numerical flux
                        // (the weak/strong form is stored in "factor")
                        solverParams.phiPsi(entity.elements[elm].edges[s], field, partialField, compField, j, factor, true, indexJ, 0,
                        					solverParams, domainDiv.nodePrec[rank]);

					}
					else // general case
					{
					    // [TO DO]: neighbours in other entities
					    // global index of the node "in front"
                       	unsigned int indexFrontJ =
                       				entity
                       					.elements[entity.elements[elm].edges[s].edgeInFront.first]
                       					.offsetInU
                       				+ entity
                       					.elements[entity.elements[elm].edges[s].edgeInFront.first]
                       					.edges[entity.elements[elm].edges[s].edgeInFront.second]
                       					.offsetInElm[entity.elements[elm].edges[s].nodeIndexEdgeInFront[j]];

                        // compute the numerical flux
                        // (the weak/strong form is stored in "factor")
                        solverParams.phiPsi(entity.elements[elm].edges[s], field, partialField, compField, j, factor, false, indexJ, 0,
                        					solverParams,  domainDiv.nodePrec[rank]);
					}
				}

				// dot product between dM and the normal
				for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
                {
                    partialField.partialIu[unk] +=
                    	entity.elements[elm].edges[s].determinantLD[0]*(
							entity.elements[elm].edges[s].normal[0]*entity.elements[elm].dM[s]*partialField.g[0][unk]
							+ entity.elements[elm].edges[s].normal[1]*entity.elements[elm].dM[s]*partialField.g[1][unk]);
                }
			}

			// add the local rhs vector to the global one
			// [TO DO]: find some Eigen function that allows to do that efficiently
			for(unsigned short unk = 0 ; unk < field.Iu.size() ; ++unk)
            {
                for(unsigned int j = 0 ; j < mesh.elementProperties.at(entity.elements[elm].elementTypeHD).nSF ; ++j)
                {
                    field.Iu[unk][entity.elements[elm].offsetInU + j - domainDiv.nodePrec[rank]] = partialField.partialIu[unk][j];
                }
            }
		}
	}
}
