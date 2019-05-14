#include <iostream>
#include <cassert>
#include <gmsh.h>
#include <mpi.h>
#include "../mpi/sendReceive.hpp"
#include "../matrices/buildMatrix.hpp"
#include "../matrices/matrix.hpp"
#include "../flux/buildFlux.hpp"
#include "../write/write.hpp"
#include "timeInteg.hpp"
#include "field.hpp"
#include "RungeKutta.hpp"
#include "../utils/utils.hpp"


/**
 * \brief Compute the increment vector of the unknown fields, for the weak form.
 * \param t Current time.
 * \param u Current solution.
 * \param field Structure that contains all the main variables.
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 */
static void Fweak(double t, Field& field, CompleteField& compField,
                    const Matrix& matrix, const Mesh& mesh,
                    const SolverParams& solverParams, const DomainDiv& domainDiv, unsigned int rank)
{
 	if(solverParams.IsSourceTerms)
        solverParams.sourceTerm(field, solverParams);

	// compute the right-hand side of the master equation (phi or psi)
 	buildFlux(mesh, field, compField, 1, t, solverParams, domainDiv, rank);

 	// compute the increment
 	for(unsigned short unk = 0 ; unk < field.DeltaU.size() ; ++unk)
    {
        field.DeltaU[unk]
            = matrix.invM*(field.Iu[unk] + matrix.Sx*field.flux[0][unk]
            				+ matrix.Sy*field.flux[1][unk]);

        if(solverParams.IsSourceTerms)
            field.DeltaU[unk]+=field.s[unk];
    }
}


/**
 * \brief Compute the increment vector of the unknown fields, for the strong form.
 * \param t Current time.
 * \param u Current solution.
 * \param field Structure that contains all the main variables.
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 */
static void Fstrong(double t, Field& field, CompleteField& compField,
                    const Matrix& matrix, const Mesh& mesh,
                    const SolverParams& solverParams, const DomainDiv& domainDiv, unsigned int rank)
{
 	if(solverParams.IsSourceTerms)
        solverParams.sourceTerm(field, solverParams);

	// compute the right-hand side of the master equation (phi or psi)
 	buildFlux(mesh, field, compField, -1, t, solverParams, domainDiv, rank);

 	// compute the increment
    for(unsigned short unk = 0 ; unk < field.DeltaU.size() ; ++unk)
    {
        field.DeltaU[unk]
            = matrix.invM*(field.Iu[unk] - matrix.Sx*field.flux[0][unk]
            				- matrix.Sy*field.flux[1][unk]);

        if(solverParams.IsSourceTerms)
            field.DeltaU[unk]+=field.s[unk];
    }
}


// see .hpp file for description
bool timeInteg(const Mesh& mesh, SolverParams& solverParams,
               const std::string& fileName, int rank, int numberProc)
{
    // Get MPI parameters
    MPI_Status status;

    if(rank == 0)
        std::cout << "Number of nodes: " << mesh.nodeData.numNodes << std::endl;

    if(numberProc >  mesh.nodeData.elementTags.size())
    {
        std::cerr << "How did you manage to end up with some processor unused ?" <<std::endl;
        return false;
    }

    DomainDiv domainDiv(numberProc);
    divideDomain(domainDiv, numberProc, mesh);

	/*******************************************************************************
	 *						            TIME STEPS 								   *
	 *******************************************************************************/
    unsigned int nTimeSteps
    	= static_cast<unsigned int>(solverParams.simTime/solverParams.timeStep);
    unsigned int nTimeStepsDtWrite
    	= static_cast<unsigned int>(solverParams.simTimeDtWrite/solverParams.timeStep);


	/*******************************************************************************
	 *						            MATRICES 								   *
	 *******************************************************************************/
	Matrix matrix;
	buildMatrix(mesh, matrix, domainDiv, rank);


	/*******************************************************************************
	 *						       WEAK AND STRONG FORM 						   *
	 *******************************************************************************/
  	//Function pointer to the used function (weak vs strong form)
  	UsedF usedF;

  	if(solverParams.solverType == "weak")
  	{
  		usedF = Fweak;
      	matrix.Sx = matrix.Sx.transpose();
      	matrix.Sy = matrix.Sy.transpose();
  	}
  	else
  	{
     	usedF = Fstrong;
  	}

  	//Initialization of the field of unknowns
  	Field field(domainDiv.node[rank], solverParams.nUnknowns, mesh.dim);
  	CompleteField compField(mesh.nodeData.numNodes, solverParams.nUnknowns, mesh.dim);
  	PartialField partialField(solverParams.nUnknowns, mesh.dim);

	/*******************************************************************************
	 *						       INITIAL CONDITION      						   *
	 *******************************************************************************/
	std::vector<double> uIC(solverParams.nUnknowns);

	for(auto entity : mesh.entities)
    {
        for(auto element : entity.elements)
        {
            for(unsigned int n = 0 ; n < element.nodeTags.size() ; ++n)
            {
                solverParams.initCondition.ibcFunc(uIC, element.nodesCoord[n], 0, field, 0, {},
					solverParams.initCondition.coefficients, solverParams.fluxCoeffs);

                for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
                    compField.u[unk](element.offsetInU + n) = uIC[unk];
            }
        }
    }

    //Each MPI Thread initialize its unknown vector from the complete one
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        for(unsigned int n = domainDiv.nodePrec[rank] ;
            n < domainDiv.nodePrec[rank] + domainDiv.node[rank] ; ++n)
        {
            field.u[unk][n - domainDiv.nodePrec[rank]]= compField.u[unk][n];
        }
    }
    //And the flux too
    solverParams.flux(field, partialField, solverParams, false);

    MPI_Barrier(MPI_COMM_WORLD);
    exchangeFlux(field, compField, domainDiv, rank, solverParams, mesh);
    MPI_Barrier(MPI_COMM_WORLD);

	/*******************************************************************************
	 *						       LAUNCH GMSH  	      						   *
	 *******************************************************************************/
    std::vector<std::string> names;
    int viewTag;
    std::string modelName;
    std::string dataType = "ElementNodeData";
    std::vector<int> elementTags = mesh.nodeData.elementTags;
    std::vector<unsigned int> elementNumNodes = mesh.nodeData.elementNumNodes;

    std::vector<std::vector<double>> uDisplay(elementNumNodes.size());
    if(rank == 0)
    {
        gmsh::initialize();
        gmsh::option::setNumber("General.Terminal", 1);
        gmsh::open(fileName);
        viewTag = gmsh::view::add("results");
        gmsh::model::list(names);
        modelName = names[0];
    }

	double t = 0.0;

	/*******************************************************************************
	 *						       INITIAL CONDITION      						   *
	 *******************************************************************************/
    if(rank == 0)
    {
		solverParams.write(uDisplay, elementNumNodes, elementTags, modelName,
                           0, 0, compField, solverParams.fluxCoeffs,
                           solverParams.whatToWrite, solverParams.viewTags);
	}

	/*******************************************************************************
	 *						       TIME INTEGRATION      						   *
	 *******************************************************************************/
	// temporary vectors (only for RK4, but I don't want to define them at each time
	// iteration)
	Field temp = field;

    //Function pointer to the used integration scheme
    IntegScheme integScheme;

    if (solverParams.timeIntType == "RK1")
    {
        integScheme = RK1;
    }
    else if (solverParams.timeIntType == "RK2")
    {
        integScheme = RK2;
    }
    else if (solverParams.timeIntType == "RK3")
    {
        integScheme = RK3;
    }
    else if (solverParams.timeIntType == "RK4")
    {
        integScheme = RK4;
    }


	// numerical integration
	unsigned int ratio, currentDecade = 0;
	for(unsigned int nbrStep = 1 ; nbrStep < nTimeSteps + 1 ;
		nbrStep++)
	{
  		// display progress
  		if(rank == 0)
        {
            ratio = int(100*double(nbrStep - 1)/double(nTimeSteps));
            if(ratio >= currentDecade)
            {
                std::cout  	<< "\r" << "Integrating: " << ratio << "%"
                            << " of the time steps done" << std::flush;
                currentDecade = ratio + 1;
            }
        }

        integScheme(t, field, partialField, compField, matrix, domainDiv, rank,
                    mesh, solverParams, temp, usedF);

        temp = field;

		// check that it does not diverge
		// assert(field.u[0].maxCoeff() <= 1E5);

		// add time step
		t += solverParams.timeStep;

        // store the results every Dt only.
		if((nbrStep % nTimeStepsDtWrite) == 0 && rank == 0)
        {
            solverParams.write(uDisplay, elementNumNodes, elementTags, modelName,
                         nbrStep, t, compField, solverParams.fluxCoeffs,
                         solverParams.whatToWrite, solverParams.viewTags);
        }
	}

	if(rank == 0)
    {
        std::cout 	<< "\r" << "Integrating: 100% of the time steps done" << std::flush
	 			<< std::endl;

		// write the results & finalize
		writeEnd(solverParams.viewTags, solverParams.whatToWrite);
		gmsh::finalize();
	}

	return true;
}
