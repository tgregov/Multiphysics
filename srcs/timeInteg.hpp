#ifndef timeInteg_hpp
#define timeInteg_hpp
#include "Mesh2D.hpp"
#include <string>


// [TO DO]: comment & describe
Eigen::VectorXd F(double t, Eigen::VectorXd& u, Eigen::VectorXd& fx, 
	Eigen::VectorXd& fy, Eigen::SparseMatrix<double> invM, 
	Eigen::SparseMatrix<double> Sx, Eigen::SparseMatrix<double> Sy, 
	unsigned int numNodes, Mesh2D& mesh, const std::string& typeForm);


bool timeInteg(Mesh2D& mesh, const std::string& scheme, const double& h, 
	const unsigned int& nbrTimeSteps, const std::string& typeForm, 
	const std::string& fileName);

#endif /* timeInteg_hpp */
