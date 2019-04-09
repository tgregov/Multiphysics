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

    field.DeltaH = invM*(IH + SxTranspose*field.FxH + SyTranspose*field.FyH);
    field.DeltauH = invM*(IuH + SxTranspose*field.FxuH + SyTranspose*field.FyuH);
    field.DeltavH = invM*(IvH + SxTranspose*field.FxvH + SyTranspose*field.FyvH);
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
	Eigen::VectorXd k1(mesh.numNodes), k2(mesh.numNodes), k3(mesh.numNodes), k4(mesh.numNodes),
	temp(mesh.numNodes);

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

		if(solverParams.timeIntType == "RK1") //(i.e. explicit Euler)
		{

			usedF(t, field, invM, Sx, Sy, mesh, solverParams.boundaryConditions);
			field.H += field.DeltaH*solverParams.timeStep;
			field.uH += field.DeltauH*solverParams.timeStep;
			field.vH += field.DeltavH*solverParams.timeStep;
		}
		/*else if(solverParams.timeIntType == "RK4")
		{
			// could be optimized
			k1 = usedF(t, u, fx, fy, invM, Sx, Sy,
						mesh, solverParams.boundaryConditions, solverParams.fluxCoeffs);

			temp = u + k1*solverParams.timeStep/2;
			k2 = usedF(t + solverParams.timeStep/2, temp, fx, fy, invM, Sx, Sy,
						mesh, solverParams.boundaryConditions, solverParams.fluxCoeffs);

			temp = u + k2*solverParams.timeStep/2;
			k3 = usedF(t + solverParams.timeStep/2, temp, fx, fy, invM, Sx, Sy,
						mesh, solverParams.boundaryConditions, solverParams.fluxCoeffs);

			temp = u + k3*solverParams.timeStep;
			k4 = usedF(t + solverParams.timeStep, temp, fx, fy, invM, Sx, Sy,
						mesh, solverParams.boundaryConditions, solverParams.fluxCoeffs);

			u += (k1 + 2*k2 + 2*k3 + k4)*solverParams.timeStep/6;

		}
		else if(solverParams.timeIntType == "RK2")
		{

			temp = u + usedF(t + solverParams.timeStep/2, u, fx, fy, invM, Sx,
							Sy, mesh, solverParams.boundaryConditions, solverParams.fluxCoeffs)
							*solverParams.timeStep/2;

			u += solverParams.timeStep * usedF(t + solverParams.timeStep/2, temp,
												fx, fy, invM, Sx, Sy, mesh,
												solverParams.boundaryConditions, solverParams.fluxCoeffs);

		}*/

		// check that it does not diverge
		// assert(field.H.maxCoeff() <= 1E5);

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

            gmsh::view::addModelData(viewTag, nbrStep, modelName,dataType, elementTags,
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
