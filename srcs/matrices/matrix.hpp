#ifndef matrix_hpp_included
#define matrix_hpp_included

#include <Eigen/Sparse>

/**
 * \struct Matrix
 * \brief A simple structure containing the matrices needed for DG-FEm.
 */
struct Matrix
{

	Eigen::SparseMatrix<double> invM;   /**< Inverse of the mass matrix */
	Eigen::SparseMatrix<double> Sx;     /**< X-stiffness matrix */
	Eigen::SparseMatrix<double> Sy;     /**< Y-stiffness matrix */
};

#endif /* matrix_hpp_included */
