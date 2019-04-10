/**
 * \file timeInteg.cpp
 * \brief Implementation of the required function to time integrate the DG-FEM equations.
 */

#include <iostream>
#include <cassert>
#include <gmsh.h>
#include "../matrices/buildM.hpp"
#include "../matrices/buildS.hpp"
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
static void Fweak(double t, Field& field, const Eigen::SparseMatrix<double>& invM,
	const Eigen::SparseMatrix<double>& SxTranspose,
	const Eigen::SparseMatrix<double>& SyTranspose,
	const Mesh& mesh,
    const std::map<std::string, ibc>& boundaries)
{
 	// compute the nodal physical fluxes
 	flux(field);

	// compute the right-hand side of the master equation (phi or psi)
	Eigen::VectorXd IH(mesh.numNodes); IH.setZero(); //[TO DO]: define this in timeInteg
	Eigen::VectorXd IuH(mesh.numNodes); IuH.setZero();
	Eigen::VectorXd IvH(mesh.numNodes); IvH.setZero();

 	buildFlux(mesh, IH, IuH, IvH, field, 1, t, boundaries);

	// compute the vector F to be integrated in time
	Eigen::VectorXd vectorF(mesh.numNodes);

    field.DeltaH = invM*(IH - SxTranspose*field.FxH - SyTranspose*field.FyH);
    field.DeltauH = invM*(IuH - SxTranspose*field.FxuH - SyTranspose*field.FyuH);
    field.DeltavH = invM*(IvH - SxTranspose*field.FxvH - SyTranspose*field.FyvH);
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
static void Fstrong(double t, Field& field, const Eigen::SparseMatrix<double>& invM,
	const Eigen::SparseMatrix<double>& Sx,
	const Eigen::SparseMatrix<double>& Sy,
	const Mesh& mesh,
    const std::map<std::string, ibc>& boundaries)
{
 	// compute the nodal physical fluxes
 	flux(field);

	// compute the right-hand side of the master equation (phi or psi)
	Eigen::VectorXd IH(mesh.numNodes); IH.setZero(); //[TO DO]: define this in timeInteg
	Eigen::VectorXd IuH(mesh.numNodes); IuH.setZero();
	Eigen::VectorXd IvH(mesh.numNodes); IvH.setZero();

 	buildFlux(mesh, IH, IuH, IvH, field, -1, t, boundaries);

	// compute the vector F to be integrated in time
	Eigen::VectorXd vectorF(mesh.numNodes);

    field.DeltaH = invM*(IH - Sx*field.FxH - Sy*field.FyH);
    field.DeltauH = invM*(IuH - Sx*field.FxuH - Sy*field.FyuH);
    field.DeltavH = invM*(IvH - Sx*field.FxvH - Sy*field.FyvH);
    //std::cout << "\n DeltauH:\n" << field.DeltauH << std::endl;
}

//Documentation in .hpp
bool timeInteg(const Mesh& mesh, const SolverParams& solverParams,
	const std::string& fileName)
{
    unsigned int nbreTimeSteps = static_cast<unsigned int>(solverParams.simTime/solverParams.timeStep);
    unsigned int nbreTimeStepsDtWrite = static_cast<unsigned int>(solverParams.simTimeDtWrite/solverParams.timeStep);

	std::vector<int> nodeTags =  getTags(mesh);

	// matrices of the DG method
  	Eigen::SparseMatrix<double> invM(mesh.numNodes, mesh.numNodes);
  	Eigen::SparseMatrix<double> Sx(mesh.numNodes, mesh.numNodes);
  	Eigen::SparseMatrix<double> Sy(mesh.numNodes, mesh.numNodes);
  	std::cout << "Building the invM matrix...";
  	buildM(mesh, invM);
  	//std::cout << "invM:\n" << invM;
  	std::cout 	<< "\rBuilding the invM matrix... 		Done" << std::flush
  				<< std::endl;
  	std::cout << "Building the Sx and Sy matrices...";
  	buildS(mesh, Sx, Sy);
  	//std::cout << "Sx:\n" << Sx;
  	//std::cout << "Sy:\n" << Sy;
  	std::cout 	<< "\rBuilding the Sx and Sy matrices... 	Done" << std::flush
  				<< std::endl;


  	//Function pointer to the used function (weak vs strong form)
  	std::function<void(double t, Field& field, const Eigen::SparseMatrix<double>& invM,
	const Eigen::SparseMatrix<double>& Sx,
	const Eigen::SparseMatrix<double>& Sy,
	const Mesh& mesh,
    const std::map<std::string, ibc>& boundaries)> usedF;

  	if(solverParams.solverType == "weak")
  	{
      	Sx = Sx.transpose();
      	Sy = Sy.transpose();
      	usedF = Fweak;
  	}
  	else
  	{
     	usedF = Fstrong;
  	}


  	Field field;
  	field.H.resize(mesh.numNodes);
  	field.uH.resize(mesh.numNodes);
  	field.vH.resize(mesh.numNodes);
  	field.FxH.resize(mesh.numNodes);
  	field.FxuH.resize(mesh.numNodes);
  	field.FxvH.resize(mesh.numNodes);
  	field.FyH.resize(mesh.numNodes);
  	field.FyuH.resize(mesh.numNodes);
  	field.FyvH.resize(mesh.numNodes);
  	field.DeltaH.resize(mesh.numNodes);
  	field.DeltauH.resize(mesh.numNodes);
  	field.DeltavH.resize(mesh.numNodes);


	//Set Initial Condition
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

	// launch gmsh
	gmsh::initialize();
	gmsh::option::setNumber("General.Terminal", 1);
	gmsh::open(fileName);
	int viewTag = gmsh::view::add("results");
	std::vector<std::string> names;
	gmsh::model::list(names);
	std::string modelName = names[0];
	std::string dataType = "ElementNodeData";

	// collect the element tags & their length
	std::vector<int> elementTags;
	std::vector<unsigned int> elementNumNodes;
	for(size_t ent = 0 ; ent < mesh.entities.size() ; ++ent)
	{
		Entity entity = mesh.entities[ent];

		for(size_t i = 0 ; i < entity.elements.size() ; ++i)
		{
			Element element = entity.elements[i];
			elementTags.push_back(element.elementTag);
			elementNumNodes.push_back(element.nodeTags.size());
		}
	}

    std::vector<std::vector<double>> uDisplay(elementNumNodes.size()); //Vector to display results

	double t = 0.0;

	//write initial condition
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


	// temporary vectors (only for RK4, but I don't want to define them at each time
	// iteration)
	Field temp = field;
	Eigen::VectorXd k1H(mesh.numNodes);
	Eigen::VectorXd k2H(mesh.numNodes);
	Eigen::VectorXd k3H(mesh.numNodes);
	Eigen::VectorXd k4H(mesh.numNodes);
	Eigen::VectorXd k1uH(mesh.numNodes);
	Eigen::VectorXd k2uH(mesh.numNodes);
	Eigen::VectorXd k3uH(mesh.numNodes);
	Eigen::VectorXd k4uH(mesh.numNodes);
	Eigen::VectorXd k1vH(mesh.numNodes);
	Eigen::VectorXd k2vH(mesh.numNodes);
	Eigen::VectorXd k3vH(mesh.numNodes);
	Eigen::VectorXd k4vH(mesh.numNodes);
	double h = solverParams.timeStep;

	// numerical integration
	unsigned int ratio, currentDecade = 0;
	for(unsigned int nbrStep = 1 ; nbrStep < nbreTimeSteps + 1 ;
		nbrStep++)
	{

  		// display progress
		ratio = int(100*double(nbrStep - 1)/double(nbreTimeSteps));
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
			usedF(t, field, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			field.H += field.DeltaH*h;
			field.uH += field.DeltauH*h;
			field.vH += field.DeltavH*h;
		}
		else if(solverParams.timeIntType == "RK2")
		{
			//Numerical time integration using the method of Runge-Kutta order 2
			usedF(t, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			k1H = temp.DeltaH*h;
			k1uH = temp.DeltauH*h;
			k1vH = temp.DeltavH*h;

			temp.H = field.H + k1H/2;
			temp.uH = field.uH + k1uH/2;
			temp.vH = field.vH + k1vH/2;

			usedF(t + h/2, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
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
			usedF(t, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			k1H = temp.DeltaH*h;
			k1uH = temp.DeltauH*h;
			k1vH = temp.DeltavH*h;

			temp.H = field.H + k1H/2;
			temp.uH = field.uH + k1uH/2;
			temp.vH = field.vH + k1vH/2;

			usedF(t + h/2, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			k2H = temp.DeltaH*h;
			k2uH = temp.DeltauH*h;
			k2vH = temp.DeltavH*h;
			
			temp.H = field.H - k1H + 2*k2H;
			temp.uH = field.H - k1uH + 2*k2uH;
			temp.vH = field.H - k1vH + 2*k2vH;

			usedF(t + h, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
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
			usedF(t, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			k1H = temp.DeltaH*h;
			k1uH = temp.DeltauH*h;
			k1vH = temp.DeltavH*h;

			temp.H = field.H + k1H/2;
			temp.uH = field.uH + k1uH/2;
			temp.vH = field.vH + k1vH/2;

			usedF(t + h/2, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			k2H = temp.DeltaH*h;
			k2uH = temp.DeltauH*h;
			k2vH = temp.DeltavH*h;
			
			temp.H = field.H + k2H/2;
			temp.uH = field.uH + k2uH/2;
			temp.vH = field.vH + k2vH/2;

			usedF(t + h/2, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			k3H = temp.DeltaH*h;
			k3uH = temp.DeltauH*h;
			k3vH = temp.DeltavH*h;

			temp.H = field.H + k3H;
			temp.uH = field.uH + k3uH;
			temp.vH = field.vH + k3vH;

			usedF(t + h, temp, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			k4H = temp.DeltaH*h;
			k4uH = temp.DeltauH*h;
			k4vH = temp.DeltavH*h;

			field.H += (k1H + 2*k2H + 2*k3H + k4H)/6;
			field.uH += (k1uH + 2*k2uH + 2*k3uH + k4uH)/6;
			field.vH += (k1vH + 2*k2vH + 2*k3vH + k4vH)/6;	
		}

		// check that it does not diverge
		//assert(field.H.maxCoeff() <= 1.5);

		// add time step
		t += solverParams.timeStep;

        // Store the results every Dt only.
		if((nbrStep % nbreTimeStepsDtWrite) == 0)
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
