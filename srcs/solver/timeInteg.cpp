#include <iostream>
#include <cassert>
#include <gmsh.h>
#include "../matrices/buildMatrix.hpp"
#include "../matrices/matrix.hpp"
#include "../flux/buildFlux.hpp"
#include "../write/write.hpp"
#include "timeInteg.hpp"
#include "field.hpp"
#include "RungeKutta.hpp"


/**
 * \brief Compute the increment vector of the unknown fields, for the weak form.
 * \param t Current time.
 * \param u Current solution.
 * \param field Structure that contains all the main variables.
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 */
static void Fweak(double t, Field& field,
                  const Matrix& matrix, const Mesh& mesh,
                  const SolverParams& solverParams)
{
 	// compute the nodal physical fluxes
 	PartialField partialField(solverParams.nUnknowns, mesh.dim);

 	solverParams.flux(field, partialField, solverParams, false);

 	if(solverParams.IsSourceTerms)
        solverParams.sourceTerm(field, solverParams);

	// compute the right-hand side of the master equation (phi or psi)
 	buildFlux(mesh, field, 1, t, solverParams);

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
static void Fstrong(double t, Field& field, const Matrix& matrix, const Mesh& mesh,
                    const SolverParams& solverParams)
{
    PartialField partialField(solverParams.nUnknowns, mesh.dim);

 	// compute the nodal physical fluxes
 	solverParams.flux(field, partialField, solverParams, false);

 	if(solverParams.IsSourceTerms)
        solverParams.sourceTerm(field, solverParams);

	// compute the right-hand side of the master equation (phi or psi)
 	buildFlux(mesh, field, -1, t, solverParams);

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
				const std::string& fileName, const std::string& resultsName)
{
        std::cout << "Number of nodes: " << mesh.nodeData.numNodes << std::endl;

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
	buildMatrix(mesh, matrix);


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
  	Field field(mesh.nodeData.numNodes, solverParams.nUnknowns, mesh.dim);

	/*******************************************************************************
	 *						       INITIAL CONDITION      						   *
	 *******************************************************************************/
	std::vector<double> uIC(solverParams.nUnknowns);

    for(auto element : mesh.elements)
    {
        for(unsigned int n = 0 ; n < element.nodeTags.size() ; ++n)
        {
            solverParams.initCondition.ibcFunc(uIC, element.nodesCoord[n], 0, field, 0, {},
                solverParams.initCondition.coefficients, solverParams.fluxCoeffs);

            for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
                field.u[unk](element.offsetInU + n) = uIC[unk];
        }
    }


	/*******************************************************************************
	 *						       LAUNCH GMSH  	      						   *
	 *******************************************************************************/
	gmsh::initialize();
	gmsh::option::setNumber("General.Terminal", 1);
	gmsh::open(fileName);
	std::vector<std::string> names;
	gmsh::model::list(names);
	std::string modelName = names[0];
	std::string dataType = "ElementNodeData";
	std::vector<int> elementTags = mesh.nodeData.elementTags;
	std::vector<unsigned int> elementNumNodes = mesh.nodeData.elementNumNodes;
    std::vector<std::vector<double>> uDisplay(elementNumNodes.size());

	double t = 0.0;

	/*******************************************************************************
	 *						       INITIAL CONDITION      						   *
	 *******************************************************************************/
    solverParams.write(uDisplay, elementNumNodes, elementTags, modelName,0, 0, field,
                         solverParams.fluxCoeffs, solverParams.whatToWrite,
                         solverParams.viewTags);


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
		ratio = int(100*double(nbrStep - 1)/double(nTimeSteps));
        if(ratio >= currentDecade)
        {
            std::cout  	<< "\r" << "Integrating: " << ratio << "%"
            			<< " of the time steps done" << std::flush;
            currentDecade = ratio + 1;
        }

        integScheme(t, field, matrix, mesh, solverParams, temp, usedF);

        temp = field;

		// check that it does not diverge
		// assert(field.u[0].maxCoeff() <= 1E5);

		// add time step
		t += solverParams.timeStep;

        // store the results every Dt only.
		if((nbrStep % nTimeStepsDtWrite) == 0)
        {
            solverParams.write(uDisplay, elementNumNodes, elementTags, modelName,
                         nbrStep, t, field, solverParams.fluxCoeffs,
                         solverParams.whatToWrite, solverParams.viewTags);
        }
	}

	std::cout 	<< "\r" << "Integrating: 100% of the time steps done" << std::flush
	 			<< std::endl;

	// write the results & finalize
    writeEnd(solverParams.viewTags, solverParams.whatToWrite, resultsName);
    gmsh::finalize();

	return true;
}
