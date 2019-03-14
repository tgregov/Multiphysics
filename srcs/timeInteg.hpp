#ifndef timeInteg_hpp_included
#define timeInteg_hpp_included

#include <string>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "Mesh2D.hpp"

// [TO DO]: comment & describe
Eigen::VectorXd F(double t, Eigen::VectorXd& u, Eigen::VectorXd& fx,
	Eigen::VectorXd& fy, Eigen::SparseMatrix<double> invM,
	Eigen::SparseMatrix<double> Sx, Eigen::SparseMatrix<double> Sy,
	unsigned int numNodes, Mesh2D& mesh, const std::string& typeForm);


bool timeInteg(Mesh2D& mesh, const std::string& scheme, const double& h,
	const unsigned int& nbrTimeSteps, const std::string& typeForm,
	const std::string& fileName);

#endif /* timeInteg_hpp */
