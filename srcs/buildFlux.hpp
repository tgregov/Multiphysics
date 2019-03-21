#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Mesh2D.hpp"
#include "bcFunction.hpp"


void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, double& C,
			const Eigen::VectorXd& u);

void flux(double& fx, double& fy, double u);


bool buildFlux(const Mesh2D& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u,
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy, const double& C,
	double factor, unsigned int numNodes, double t, const std::map<std::string, bc>& boundaries);

#endif /* buildFlux_hpp */
