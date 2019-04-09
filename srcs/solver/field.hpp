#ifndef field_hpp_included
#define field_hpp_included

#include <Eigen/Dense>

struct Field
{
	Eigen::VectorXd H;
	Eigen::VectorXd uH;
	Eigen::VectorXd vH;

	Eigen::VectorXd FxH, 	FyH;
	Eigen::VectorXd FxuH, 	FyuH; 
	Eigen::VectorXd FxvH, 	FyvH;
};
