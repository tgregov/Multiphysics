#include <iostream>
#include <string>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <gmsh.h>
#include "Mesh2D.hpp"
#include "buildM.hpp"
#include "buildS.hpp"
#include "buildFlux.hpp"


Eigen::VectorXd F(double t, Eigen::VectorXd& u, Eigen::VectorXd& fx, 
	Eigen::VectorXd& fy, Eigen::SparseMatrix<double> invM, 
	Eigen::SparseMatrix<double> Sx, Eigen::SparseMatrix<double> Sy, 
	unsigned int numNodes, Mesh2D& mesh, const std::string& typeForm)
{
	
 	// compute the nodal physical fluxes
 	double C;
 	flux(fx, fy, C, u);

	// compute the right-hand side of the master equation (phi or psi)
	Eigen::VectorXd I(numNodes); I.setZero(); //[TO DO]: define this in timeInteg
 	buildFlux(mesh, I, u, fx, fy, C, typeForm, numNodes);

	// compute the vector F to be integrated in time
	Eigen::VectorXd vectorF(numNodes);
	if(!typeForm.compare("strong"))
	{
		vectorF = invM*(I - Sx*fx - Sy*fy);
	} 
	else if(!typeForm.compare("weak")) 
	{
		vectorF = invM*(I + Sx.transpose()*fx + Sy.transpose()*fy);
	} 
	else 
	{
		std::cerr 	<< "The form  " << typeForm  << "does not exist !" << std::endl;
	}

	return vectorF;
}


bool timeInteg(Mesh2D& mesh, const std::string& scheme, const double& h, 
	const unsigned int& nbrTimeSteps, const std::string& typeForm, 
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

	// display the results
	std::cout << "Matrix [M]:\n" << M << std::endl;
	std::cout << "Matrix [M^-1]:\n" << invM << std::endl;
    std::cout << "Matrix [Sx]:\n" << Sx << std::endl;
    std::cout << "Matrix [Sy]:\n" << Sy << std::endl;

    // initial condition [TO DO]: generalize for u != 0
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
    std::string dataType = "NodeData";
	
	// numerical integration
	if(!scheme.compare("RK1")) // Runge-Kutta of order 1 (i.e. explicit Euler)
	{ 
		double t = 0.0;
		
		std::vector<std::vector<double>> uDisplay(numNodes);

		for(unsigned int count = 0; count < u.size(); ++count)
		{
			uDisplay[count].resize(1);
			uDisplay[count][0] = u[count];
		}
		gmsh::view::addModelData(viewTag, 0, modelName, dataType, nodeTags,
			uDisplay, t);

		for(unsigned int nbrStep = 1 ; nbrStep < nbrTimeSteps + 1 ; nbrStep++)
		{
			
			std::cout << "[Time step: " << nbrStep << "]" << std::endl;

			u += F(t, u, fx, fy, invM, Sx, Sy, numNodes, mesh, typeForm)*h;
			t += h;
			

			for(unsigned int count = 0; count < u.size(); ++count)
			{
				uDisplay[count].resize(1);
				uDisplay[count][0] = u[count];
			}
			gmsh::view::addModelData(viewTag, nbrStep, modelName, dataType, nodeTags,
				uDisplay, t);
		}
	} else if(!scheme.compare("RK4")){ // Runge-Kutta of order 4
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
	} else{
		std::cerr 	<< "The integration scheme " << scheme 
					<< " is not implemented" << std::endl;
		return false;
	}
	
	// write the results & finalize
    gmsh::view::write(viewTag, std::string("results.msh"));
    gmsh::finalize();

	return true;
}