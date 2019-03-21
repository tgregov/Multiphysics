#ifndef timeInteg_hpp_included
#define timeInteg_hpp_included

#include <string>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "./mesh/Mesh2D.hpp"

// [TO DO]: comment & describe
Eigen::VectorXd Fweak(double t, Eigen::VectorXd& u, Eigen::VectorXd& fx,
	Eigen::VectorXd& fy, const Eigen::SparseMatrix<double>& invM,
	const Eigen::SparseMatrix<double>& SxTranspose, 
	const Eigen::SparseMatrix<double>& SyTranspose,
	unsigned int numNodes, const Mesh2D& mesh);

Eigen::VectorXd Fstrong(double t, Eigen::VectorXd& u, Eigen::VectorXd& fx,
	Eigen::VectorXd& fy, const Eigen::SparseMatrix<double>& invM,
	const Eigen::SparseMatrix<double>& SxTranspose, 
	const Eigen::SparseMatrix<double>& SyTranspose,
	unsigned int numNodes, const Mesh2D& mesh);

bool timeInteg(Mesh2D& mesh, const std::string& scheme, const double& h,
	const unsigned int& nbrTimeSteps, const std::string& typeForm,
	const std::string& fileName);

#endif /* timeInteg_hpp */
