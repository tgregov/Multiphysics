#include <iostream>
#include <cassert>
#include <gmsh.h>
#include <string>
#include "../matrices/buildMatrix.hpp"
#include "../matrices/matrix.hpp"
#include "../flux/buildFlux.hpp"
#include "timeInteg.hpp"
#include "field.hpp"
#include "rungeKutta.hpp"


//typedef to lighten the notations
typedef std::function<void(double, Field&, PartialField&, const Matrix&,
         const Mesh&, const SolverParams&, Field&, UsedF)> IntegScheme;



/**
 * \brief Compute the increment vector of the unknown fields, for the weak form.
 * \param t Current time.
 * \param u Current solution.
 * \param field Field that contains all the main variables.
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 */
static void Fweak(double t, Field& field, PartialField& partialField, const Matrix& matrix, const Mesh& mesh,
                  const SolverParams& solverParams)
{
 	// compute the nodal physical fluxes
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
 * \param field Field that contains all the main variables.
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 */
static void Fstrong(double t, Field& field, PartialField& partialField, const Matrix& matrix, const Mesh& mesh,
                    const SolverParams& solverParams)
{
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
bool timeInteg(const Mesh& mesh, const SolverParams& solverParams,
				const std::string& fileName, const std::string& resultsName)
{

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
                    field.u[unk](element.offsetInU + n) = uIC[unk];
            }
        }
    }



	/*******************************************************************************
	 *						       LAUNCH GMSH  	      						   *
	 *******************************************************************************/
	gmsh::initialize();
	gmsh::option::setNumber("General.Terminal", 1);
	gmsh::open(fileName);
	int viewTag1 = gmsh::view::add("Height");
    int viewTag2 = gmsh::view::add("Velocity Field");
	std::vector<std::string> names;
	gmsh::model::list(names);
	std::string modelName = names[0];
	std::string dataType = "ElementNodeData";
	std::vector<int> elementTags = mesh.nodeData.elementTags;
	std::vector<unsigned int> elementNumNodes = mesh.nodeData.elementNumNodes;
    std::vector<std::vector<double>> uDisplay1(elementNumNodes.size());
    std::vector<std::vector<double>> uDisplay2(elementNumNodes.size());

	double t = 0.0;

	/*******************************************************************************
	 *						       INITIAL CONDITION      						   *
	 *******************************************************************************/
	std::cout << "coucou1" << std::endl;
    unsigned int index = 0;

	for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
	{
		std::vector<double> tempWrite1(elementNumNodes[count]);
        std::vector<double> tempWrite2(3*elementNumNodes[count]);
		for(unsigned int node = 0 ; node < elementNumNodes[count] ; ++node)
		{
			//Height
            tempWrite1[node]=field.u[0][index];
            
            //Speed 
            // tempWrite2[3*node]=field.u[1][index]/field.u[0][index];
            // tempWrite2[3*node+1]=field.u[2][index]/field.u[0][index];
            // tempWrite2[3*node+2]=0;

			++index;
		}

		uDisplay1[count] = std::move(tempWrite1);
        //uDisplay2[count] = std::move(tempWrite2);
	}

	gmsh::view::addModelData(viewTag1, 0, modelName, dataType, elementTags,
	                         uDisplay1, t, 1);
    //gmsh::view::addModelData(viewTag2, 0, modelName, dataType, elementTags,
    //                         uDisplay2, t, 3);

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

        //Call to RK
        integScheme(t, field, partialField, matrix, mesh, solverParams, temp, usedF);
        
        //Update of temporary field
        temp = field;

		// check that it does not diverge
		// assert(field.u[0].maxCoeff() <= 1E5);

		// add time step
		t += solverParams.timeStep;

        // store the results every Dt only.
		if((nbrStep % nTimeStepsDtWrite) == 0)
        {
            unsigned int offset = 0;
            for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
            {
                std::vector<double> tempWrite1(elementNumNodes[count]);
                std::vector<double> tempWrite2(3*elementNumNodes[count]);
                for (unsigned int countLocal = 0; countLocal < elementNumNodes[count];
                    ++countLocal)
                {
                    //Height
                    tempWrite1[countLocal] = field.u[0][countLocal+offset];

                    //Velocity (vector of 2 components)
                    // tempWrite2[3*countLocal] = field.u[1][countLocal+offset]/field.u[0][countLocal+offset];
                    // tempWrite2[3*countLocal+1] = field.u[2][countLocal+offset]/field.u[0][countLocal+offset];
                    // tempWrite2[3*countLocal+2] = 0;
                }

                offset += elementNumNodes[count];
                uDisplay1[count] = std::move(tempWrite1);
                //uDisplay2[count] = std::move(tempWrite2);
            }

            gmsh::view::addModelData(viewTag1, nbrStep, modelName, dataType, elementTags,
                uDisplay1, t, 1);

            // gmsh::view::addModelData(viewTag2, nbrStep, modelName, dataType, elementTags,
            //     uDisplay2, t, 3);
        }
	}

	std::cout 	<< "\r" << "Integrating: 100% of the time steps done" << std::flush
	 			<< std::endl;


    gmsh::view::write(viewTag1, resultsName);
    //gmsh::view::write(viewTag2, "./results/Velocity.msh");
    gmsh::finalize();


	return true;
}
