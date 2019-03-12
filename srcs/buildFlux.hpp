#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <Eigen/Dense>
#include "readMesh.hpp"

void flux(Eigen::VectorXd& fx, Eigen::VectorXd& fy, double& C,
			const Eigen::VectorXd& u);
bool buildFlux(Mesh2D& mesh, Eigen::VectorXd& I,
				const Eigen::VectorXd& u, const std::string& typeForm, unsigned int numNodes);

#endif /* buildFlux_hpp */
