#ifndef buildM_hpp
#define buildM_hpp

#include <Eigen/Sparse>
#include "Mesh2D.hpp"


/**
 * \brief Function that builds the [M] matrix of the DG method. This matrix is
 * defined as M_ij = sum_k{w_k*l_i(x_k)*l_j(x_k)*det[J](x_k)}, where the sum is done
 * over the Gauss points (GP). Since the components w_k*l_i(x_k)*l_j(x_k) are already
 * computed in the mesh, it suffices to get the determinant of the elements, 
 * calculate this sum, and store the result in a sparse matrix
 * \param mesh2D The structure that contains the mesh.
 * \param M The Eigen::SparseMatrix in which the matrix components will be stored.
 */
void buildM(Mesh2D& mesh2D, Eigen::SparseMatrix<double>& M);

#endif /* buildM_hpp */
