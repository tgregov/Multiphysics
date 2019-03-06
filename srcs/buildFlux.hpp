#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <Eigen/Dense>
#include "readMesh.hpp"

void flux(Eigen::VectorXd<double>& fx, Eigen::VectorXd<double>& fy, double& C, 
			const Eigen::VectorXd<double>& u);
bool buildFlux(const MeshParams& meshParams, Eigen::VectorXd<double>& I,
				const Eigen::VectorXd<double>& u, const std::string& typeForm);

#endif /* buildFlux_hpp */
