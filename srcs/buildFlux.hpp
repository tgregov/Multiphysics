#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <Eigen/Dense>


void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, double& C,
			const Eigen::VectorXd& u);

void flux(double& fx, double& fy, double u);


bool buildFlux(const Mesh2D& mesh, Eigen::VectorXd& I, const Eigen::VectorXd& u,
	const Eigen::VectorXd& fx, const Eigen::VectorXd& fy, const double& C,
	const std::string& typeForm, unsigned int numNodes, double t);

#endif /* buildFlux_hpp */
