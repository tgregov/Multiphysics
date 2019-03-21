/**
 * \file timeInteg.cpp
 * \brief Implementation of the required function to time integrate the DG-FEM equations.
 */

#include <iostream>
#include <gmsh.h>
#include "./matrices/buildM.hpp"
#include "./matrices/buildS.hpp"
#include "buildFlux.hpp"
#include "timeInteg.hpp"
#include "bcFunction.hpp"

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
static Eigen::VectorXd Fweak(double t, Eigen::VectorXd& u, Eigen::VectorXd& fx,
	Eigen::VectorXd& fy, const Eigen::SparseMatrix<double>& invM,
	const Eigen::SparseMatrix<double>& SxTranspose, 
	const Eigen::SparseMatrix<double>& SyTranspose,
	unsigned int numNodes, const Mesh2D& mesh,
  const std::map<std::string, bc>& boundaries)
{

 	// compute the nodal physical fluxes
 	double C;
 	flux(fx, fy, C, u);

	// compute the right-hand side of the master equation (phi or psi)
	Eigen::VectorXd I(numNodes); I.setZero(); //[TO DO]: define this in timeInteg

 	buildFlux(mesh, I, u, fx, fy, C, 1, numNodes, t, boundaries);

	// compute the vector F to be integrated in time
	Eigen::VectorXd vectorF(numNodes);

    vectorF = invM*(I + SxTranspose*fx + SyTranspose*fy);

	return vectorF;
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
static Eigen::VectorXd Fstrong(double t, Eigen::VectorXd& u, Eigen::VectorXd& fx,
	Eigen::VectorXd& fy, const Eigen::SparseMatrix<double>& invM,
	const Eigen::SparseMatrix<double>& Sx, const Eigen::SparseMatrix<double>& Sy,
	unsigned int numNodes, const Mesh2D& mesh, const std::map<std::string, bc>& boundaries)
{

 	// compute the nodal physical fluxes
 	double C;
 	flux(fx, fy, C, u);

	// compute the right-hand side of the master equation (phi or psi)
	Eigen::VectorXd I(numNodes); I.setZero(); //[TO DO]: define this in timeInteg

 	buildFlux(mesh, I, u, fx, fy, C, -1, numNodes, t, boundaries);

	// compute the vector F to be integrated in time
	Eigen::VectorXd vectorF(numNodes);

    vectorF = invM*(I - Sx*fx - Sy*fy);

	return vectorF;
}

//Documentation in .hpp
bool timeInteg(const Mesh2D& mesh, const SolverParams& solverParams,
	const std::string& fileName)
{

	// number of nodes and tags of the problem
	unsigned int numNodes = getNumNodes(mesh);
	std::vector<int> nodeTags =  getTags(mesh);

	// matrices of the DG method
  Eigen::SparseMatrix<double> M(numNodes, numNodes);
  Eigen::SparseMatrix<double> Sx(numNodes, numNodes);
  Eigen::SparseMatrix<double> Sy(numNodes, numNodes);
  buildM(mesh, M);
  buildS(mesh, Sx, Sy);

  //Function pointer to the used function (weak vs strong form)
  std::function<Eigen::VectorXd(double t, Eigen::VectorXd& u,
                                Eigen::VectorXd& fx, Eigen::VectorXd& fy,
                                const Eigen::SparseMatrix<double>& invM,
                                const Eigen::SparseMatrix<double>& SxTranspose,
                                const Eigen::SparseMatrix<double>& SyTranspose,
                                unsigned int numNodes, const Mesh2D& mesh,
                                const std::map<std::string, bc>& boundaries)> usedF;

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

	// invert [M]
	Eigen::SparseMatrix<double> eye(numNodes, numNodes);
	eye.setIdentity();
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solverM;
	solverM.compute(M);
	if(solverM.info() != Eigen::Success)
	{
		return false;
	}
	Eigen::SparseMatrix<double> invM(numNodes, numNodes);
	invM = solverM.solve(eye);

  // initial condition [TO DO]: use param.dat and bc struct
  Eigen::VectorXd u(numNodes); u.setZero();

  // vectors of physical flux
  Eigen::VectorXd fx(numNodes);
  Eigen::VectorXd fy(numNodes);

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
      Entity2D entity = mesh.entities[ent];

      for(size_t i = 0 ; i < entity.elements.size() ; ++i)
      {
        Element2D element = entity.elements[i];
        elementTags.push_back(element.elementTag);
        // [TO DO]: only works for T3 :(
        elementNumNodes.push_back(element.edges.size());
      }
  }

  double t = 0.0;

  //write initial condition
  std::vector<std::vector<double>> uDisplay;
  unsigned int index = 0;

  for(size_t count = 0 ; count < elementTags.size() ; ++count)
  {

      std::vector<double> temp;
      for(unsigned int node = 0 ; node < elementNumNodes[count] ; ++node)
      {
          temp.push_back(u[index]);
          ++index;
      }

      uDisplay.push_back(temp);
  }

  gmsh::view::addModelData(viewTag, 0, modelName, dataType, elementTags,
      uDisplay, t, 1);

	// numerical integration
	if(solverParams.timeIntType == "RK1") // Runge-Kutta of order 1 (i.e. explicit Euler)
	{
		for(unsigned int nbrStep = 1 ; nbrStep < solverParams.nbrTimeSteps + 1 ; nbrStep++)
		{
			std::cout << "[Time step: " << nbrStep << "]" << std::endl;

			u += usedF(t, u, fx, fy, invM, Sx, Sy, numNodes, mesh, solverParams.boundaryConditions)*solverParams.timeStep;

			/*
			for(size_t i = 0; i < u.size(); i++){
				std::cout << "u after TS [" << i << "]: " <<  u[i] << std::endl;
			}
			*/

			t += solverParams.timeStep;

			std::vector<std::vector<double>> uDisplay;

			for(unsigned int count = 0; count < u.size()/3; ++count)
			{
				std::vector<double> temp;
				temp.push_back(u[3*count]);
				temp.push_back(u[3*count+1]);
				temp.push_back(u[3*count+2]);
				uDisplay.push_back(temp);
			}

			gmsh::view::addModelData(viewTag, nbrStep, modelName,dataType, 
				elementTags, uDisplay, t, 1);
		}
	}
	else if(solverParams.timeIntType == "RK4"){ // Runge-Kutta of order 4
		/*

		// k1, k2, k3, k4 are the same dimensions as u[i] => get the best type
		// in Eigen
		Eigen::VectorXd k1, k2, k3, k4;

		for(int i = 1; i < nbrTimeSteps+1; i++)
		{
			// To change to use the abilities of Eigen
			k1 = F(t, u[i-1], M, S, uBC);
			k2 = F(t + h/2, u[i-1] + h*k1/2, invM, Sx, Sy);
			k3 = F(t + h/2, u[i-1] + h*k2/2, M, S, uBC);
			k4 = F(t + h, u[i-1] + h*k3, M, S, uBC);

			u[i] = u[i-1] + (k1 + 2*k2 + 2*k3 + k4)*h/6;
			t += h;
		}*/
	}

	// write the results & finalize
    gmsh::view::write(viewTag, std::string("results.msh"));
    gmsh::finalize();

	return true;
}
