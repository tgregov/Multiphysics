#ifndef field_hpp_included
#define field_hpp_included

#include <Eigen/Dense>

struct Field
{
	// solution fields
	Eigen::VectorXd H;
	Eigen::VectorXd uH;
	Eigen::VectorXd vH;

	// physical flux fields
	Eigen::VectorXd FxH, 	FyH;
	Eigen::VectorXd FxuH, 	FyuH; 
	Eigen::VectorXd FxvH, 	FyvH;

	// time-integration increment
	Eigen::VectorXd DeltaH;
	Eigen::VectorXd DeltauH;
	Eigen::VectorXd DeltavH;

	// rhs fields
	Eigen::VectorXd IH;
	Eigen::VectorXd IuH;
	Eigen::VectorXd IvH;

};
#endif