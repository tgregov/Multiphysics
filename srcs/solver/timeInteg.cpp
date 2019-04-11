/**
 * \file timeInteg.cpp
 * \brief Implementation of the required function to time integrate the DG-FEM equations.
 */

#include <iostream>
#include <cassert>
#include <gmsh.h>
#include "../matrices/buildMatrix.hpp"
#include "../matrices/matrix.hpp"
#include "buildFlux.hpp"
#include "timeInteg.hpp"
#include "../params/ibvFunction.hpp"
#include "field.hpp"

/**
 * \brief Compute the Du vector.
 * \param t Current time.
 * \param u Current solution.
 * \param fx Vector of flux along x.
 * \param fy Vector of flux along y.
 * \param invM Inverse of the mass matrix.
 * \param SxTranspose Transpose of the x stiffness matrix.
 * \param SyTranspose Transpose of the y stiffness matrix.
 * \param numNodes Number of nodes in the mesh.
 * \param mesh Mesh representing the domain.
 * \param boundaries Map to access mathematical function
 * \return DU vector.
 * for each BC at each boundaries.
 */
static void Fweak(double t, Field& field, const Matrix& matrix, const Mesh& mesh,
    				const std::map<std::string, ibc>& boundaries)
{
 	// compute the nodal physical fluxes
 	flux(field);

	// compute the right-hand side of the master equation (phi or psi)
 	buildFlux(mesh, field, 1, t, boundaries);

 	// compute the increment
    field.DeltaH 
    	= matrix.invM*(field.IH + matrix.Sx*field.FxH + matrix.Sy*field.FyH);
    field.DeltauH 
    	= matrix.invM*(field.IuH + matrix.Sx*field.FxuH + matrix.Sy*field.FyuH);
    field.DeltavH 
    	= matrix.invM*(field.IvH + matrix.Sx*field.FxvH + matrix.Sy*field.FyvH);
}

/**
 * \brief Compute the Du vector.
 * \param t Current time.
 * \param u Current solution.
 * \param fx Vector of flux along x.
 * \param fy Vector of flux along y.
 * \param invM Inverse of the mass matrix.
 * \param Sx x stiffness matrix.
 * \param Sy y stiffness matrix.
 * \param numNodes Number of nodes in the mesh.
 * \param mesh Mesh representing the domain.
 * \param boundaries Map to access mathematical function
 * \return DU vector.
 * for each BC at each boundaries.
 */
static void Fstrong(double t, Field& field, const Matrix& matrix, const Mesh& mesh,
    				const std::map<std::string, ibc>& boundaries)
{
 	// compute the nodal physical fluxes
 	flux(field);

	// compute the right-hand side of the master equation (phi or psi)
 	buildFlux(mesh, field, -1, t, boundaries);

 	// compute the increment
    field.DeltaH 
    	= matrix.invM*(field.IH - matrix.Sx*field.FxH - matrix.Sy*field.FyH);
    field.DeltauH 
    	= matrix.invM*(field.IuH - matrix.Sx*field.FxuH - matrix.Sy*field.FyuH);
    field.DeltavH 
    	= matrix.invM*(field.IvH - matrix.Sx*field.FxvH - matrix.Sy*field.FyvH);
}

//Documentation in .hpp
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
  	std::function<void(double t, Field& field, const Matrix& matrix,
	const Mesh& mesh,
    const std::map<std::string, ibc>& boundaries)> usedF;

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
  	Field field(mesh.nodeData.numNodes);


	/*******************************************************************************
	 *						       INITIAL CONDITION      						   *
	 *******************************************************************************/
	field.uH.setZero();
	field.vH.setZero();

	for(auto entity : mesh.entities)
    {
        for(auto element : entity.elements)
        {
            for(unsigned int n = 0 ; n < element.nodeTags.size() ; ++n)
            {
                field.H(element.offsetInU + n) =
                solverParams.initCondition.ibcFunc(element.nodesCoord[n], 0, 0,
					solverParams.initCondition.coefficients);
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
			temp[node]=field.H[index];
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
	Eigen::VectorXd k1H(mesh.nodeData.numNodes);
	Eigen::VectorXd k2H(mesh.nodeData.numNodes);
	Eigen::VectorXd k3H(mesh.nodeData.numNodes);
	Eigen::VectorXd k4H(mesh.nodeData.numNodes);
	Eigen::VectorXd k1uH(mesh.nodeData.numNodes);
	Eigen::VectorXd k2uH(mesh.nodeData.numNodes);
	Eigen::VectorXd k3uH(mesh.nodeData.numNodes);
	Eigen::VectorXd k4uH(mesh.nodeData.numNodes);
	Eigen::VectorXd k1vH(mesh.nodeData.numNodes);
	Eigen::VectorXd k2vH(mesh.nodeData.numNodes);
	Eigen::VectorXd k3vH(mesh.nodeData.numNodes);
	Eigen::VectorXd k4vH(mesh.nodeData.numNodes);
	double h = solverParams.timeStep;

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

		if(solverParams.timeIntType == "RK1")
		{
			//Numerical time integration using the method of Runge-Kutta order 1
			// i.e. explicit Euler
			usedF(t, field, matrix, mesh, solverParams.boundaryConditions);
			field.H += field.DeltaH*h;
			field.uH += field.DeltauH*h;
			field.vH += field.DeltavH*h;
		}
		else if(solverParams.timeIntType == "RK2")
		{
			//Numerical time integration using the method of Runge-Kutta order 2
			usedF(t, temp, matrix, mesh, solverParams.boundaryConditions);
			k1H = temp.DeltaH*h;
			k1uH = temp.DeltauH*h;
			k1vH = temp.DeltavH*h;

			temp.H = field.H + k1H/2;
			temp.uH = field.uH + k1uH/2;
			temp.vH = field.vH + k1vH/2;

			usedF(t + h/2, temp, matrix, mesh, solverParams.boundaryConditions);
			k2H = temp.DeltaH*h;
			k2uH = temp.DeltauH*h;
			k2vH = temp.DeltavH*h;

			field.H += k2H;
			field.uH += k2uH;
			field.vH += k2vH;
		}
		else if(solverParams.timeIntType == "RK3")
		{
			//Numerical time integration using the method of Runge-Kutta order 3
			usedF(t, temp, matrix, mesh, solverParams.boundaryConditions);
			k1H = temp.DeltaH*h;
			k1uH = temp.DeltauH*h;
			k1vH = temp.DeltavH*h;

			temp.H = field.H + k1H/2;
			temp.uH = field.uH + k1uH/2;
			temp.vH = field.vH + k1vH/2;

			usedF(t + h/2, temp, matrix, mesh, solverParams.boundaryConditions);
			k2H = temp.DeltaH*h;
			k2uH = temp.DeltauH*h;
			k2vH = temp.DeltavH*h;
			
			temp.H = field.H - k1H + 2*k2H;
			temp.uH = field.uH - k1uH + 2*k2uH;
			temp.vH = field.vH - k1vH + 2*k2vH;

			usedF(t + h, temp, matrix, mesh, solverParams.boundaryConditions);
			k3H = temp.DeltaH*h;
			k3uH = temp.DeltauH*h;
			k3vH = temp.DeltavH*h;

			field.H += (k1H + 4*k2H + k3H)/6;
			field.uH += (k1uH + 4*k2uH + k3uH)/6;
			field.vH += (k1vH + 4*k2vH + k3vH)/6;
		}
		else if(solverParams.timeIntType == "RK4")
		{
			//Numerical time integration using the method of Runge-Kutta order 4
			usedF(t, temp, matrix, mesh, solverParams.boundaryConditions);
			k1H = temp.DeltaH*h;
			k1uH = temp.DeltauH*h;
			k1vH = temp.DeltavH*h;

			temp.H = field.H + k1H/2;
			temp.uH = field.uH + k1uH/2;
			temp.vH = field.vH + k1vH/2;

			usedF(t + h/2, temp, matrix, mesh, solverParams.boundaryConditions);
			k2H = temp.DeltaH*h;
			k2uH = temp.DeltauH*h;
			k2vH = temp.DeltavH*h;
			
			temp.H = field.H + k2H/2;
			temp.uH = field.uH + k2uH/2;
			temp.vH = field.vH + k2vH/2;

			usedF(t + h/2, temp, matrix, mesh, solverParams.boundaryConditions);
			k3H = temp.DeltaH*h;
			k3uH = temp.DeltauH*h;
			k3vH = temp.DeltavH*h;

			temp.H = field.H + k3H;
			temp.uH = field.uH + k3uH;
			temp.vH = field.vH + k3vH;

			usedF(t + h, temp, matrix, mesh, solverParams.boundaryConditions);
			k4H = temp.DeltaH*h;
			k4uH = temp.DeltauH*h;
			k4vH = temp.DeltavH*h;

			field.H += (k1H + 2*k2H + 2*k3H + k4H)/6;
			field.uH += (k1uH + 2*k2uH + 2*k3uH + k4uH)/6;
			field.vH += (k1vH + 2*k2vH + 2*k3vH + k4vH)/6;	
		}

		// check that it does not diverge
		// assert(field.H.maxCoeff() <= 1E5);

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
                    temp[countLocal] = field.H[countLocal+offset];
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
