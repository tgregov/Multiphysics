#include <iostream>
#include <cassert>
#include <gmsh.h>
#include "../matrices/buildMatrix.hpp"
#include "../matrices/matrix.hpp"
#include "../flux/buildFlux.hpp"
#include "timeInteg.hpp"
#include "field.hpp"

//typedef to lighten the notations
typedef std::function<void (double, Field &, PartialField &, const Matrix &,
      const Mesh &, const SolverParams &)> UsedF;

typedef std::function<void(double, Field&, PartialField&, const Matrix&,
         const Mesh&, const SolverParams&, Field&, UsedF)> IntegScheme;


/**
 * \brief Compute the numerical time integration using the method of Runge-Kutta order 1
    i.e. explicit Euler
 * \param t Current time.
 * \param field Field that contains all the main variables.
 * \param partialField Field that contains all the main variables (private)
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 * \param temp temporary Field needed to compute the different k's
 * \param usedF pointer to the function Fweak or Fstrong
 */
static void RK1(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{
    usedF(t, field, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
        field.u[unk] += field.DeltaU[unk]*solverParams.timeStep;
}


/**
 * \brief Compute the numerical time integration using the method of Runge-Kutta order 2
    i.e. explicit Euler
 * \param t Current time.
 * \param field Field that contains all the main variables.
 * \param partialField Field that contains all the main variables (private)
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 * \param temp temporary Field needed to compute the different k's
 * \param usedF pointer to the function Fweak or Fstrong
 */
static void RK2(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{

    double h = solverParams.timeStep;

    usedF(t, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += field.k2[unk];
    }
}

/**
 * \brief Compute the numerical time integration using the method of Runge-Kutta order 3
    i.e. explicit Euler
 * \param t Current time.
 * \param field Field that contains all the main variables.
 * \param partialField Field that contains all the main variables (private)
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 * \param temp temporary Field needed to compute the different k's
 * \param usedF pointer to the function Fweak or Fstrong
 */
static void RK3(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{
    double h = solverParams.timeStep;

    usedF(t, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] - field.k1[unk] + 2*field.k2[unk];
    }

    usedF(t + h, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < 3 ; ++unk)
    {
        field.k3[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += (field.k1[unk] + 4*field.k2[unk] + field.k3[unk])/6;
    }
}

/**
 * \brief Compute the numerical time integration using the method of Runge-Kutta order 4
    i.e. explicit Euler
 * \param t Current time.
 * \param field Field that contains all the main variables.
 * \param partialField Field that contains all the main variables (private)
 * \param matrix Structure that contains the matrices of the DG method.
 * \param mesh Mesh representing the domain.
 * \param solverParams Parameters of the solver.
 * \param temp temporary Field needed to compute the different k's
 * \param usedF pointer to the function Fweak or Fstrong
 */
static void RK4(double t, Field& field, PartialField& partialField, const Matrix& matrix,
         const Mesh& mesh, const SolverParams& solverParams, Field& temp, UsedF usedF)
{
    double h = solverParams.timeStep;

    usedF(t, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k1[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k1[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k2[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k2[unk]/2;
    }

    usedF(t + h/2, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k3[unk] = temp.DeltaU[unk]*h;
        temp.u[unk] = field.u[unk] + field.k3[unk];
    }

    usedF(t + h, temp, partialField, matrix, mesh, solverParams);
    for(unsigned short unk = 0 ; unk < solverParams.nUnknowns ; ++unk)
    {
        field.k4[unk] = temp.DeltaU[unk]*h;
        field.u[unk] += (field.k1[unk] + 2*field.k2[unk] + 2*field.k3[unk] + field.k4[unk])/6;
    }
}

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
				const std::string& fileName)
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
  	std::function<void(double t, Field& field, PartialField& partialField, const Matrix& matrix,
	const Mesh& mesh,
    const SolverParams& solverParams)> usedF;

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
	int viewTag = gmsh::view::add("results");
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
	unsigned int index = 0;

	for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
	{
		std::vector<double> temp(elementNumNodes[count]);
		for(unsigned int node = 0 ; node < elementNumNodes[count] ; ++node)
		{
			temp[node]=field.u[0][index];
			++index;
		}

		uDisplay[count] = std::move(temp);
	}

	gmsh::view::addModelData(viewTag, 0, modelName, dataType, elementTags,
	                         uDisplay, t, 1);


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


        integScheme(t, field, partialField, matrix, mesh, solverParams, temp, usedF);


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
                std::vector<double> temp(elementNumNodes[count]);
                for (unsigned int countLocal = 0; countLocal < elementNumNodes[count];
                    ++countLocal)
                {
                    temp[countLocal] = field.u[0][countLocal+offset];
                }
                offset += elementNumNodes[count];
                uDisplay[count] = std::move(temp);
            }

            gmsh::view::addModelData(viewTag, nbrStep, modelName, dataType, elementTags,
                uDisplay, t, 1);
        }
	}

	std::cout 	<< "\r" << "Integrating: 100% of the time steps done" << std::flush
	 			<< std::endl;

	// write the results & finalize
    gmsh::view::write(viewTag, std::string("results.msh"));
    gmsh::finalize();


	return true;
}
