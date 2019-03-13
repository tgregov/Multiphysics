#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <Eigen/Dense>


void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, double& C,
			const Eigen::VectorXd& u);


bool buildFlux(Mesh2D& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u, 
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy, const double& C,
	const std::string& typeForm, unsigned int numNodes);

#endif /* buildFlux_hpp */
