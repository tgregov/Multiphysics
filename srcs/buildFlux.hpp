#ifndef buildFlux_hpp
#define buildFlux_hpp

#include <Eigen/Dense>
#include "readMesh.hpp"

void flux(Eigen::VectorXd<double>& fx, Eigen::VectorXd<double>& fy, double& C, 
			const Eigen::VectorXd<double>& u);
void buildFlux(MeshParams& meshParams, Eigen::VectorXd<double>& I);

#endif /* buildFlux_hpp */
