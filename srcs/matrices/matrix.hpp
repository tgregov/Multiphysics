#ifndef matrix_hpp_included
#define matrix_hpp_included

#include <Eigen/Sparse>

struct Matrix
{

	Eigen::SparseMatrix<double> invM;
	Eigen::SparseMatrix<double> Sx;
	Eigen::SparseMatrix<double> Sy;
};

#endif