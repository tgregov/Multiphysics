#ifndef buildS_hpp
#define buildS_hpp

#include <Eigen/Sparse>
#include "../mesh/Mesh.hpp"


/**
 * \brief Function that builds the [Sx], [Sy] matrix of the DG method. This matrix is
 * defined as
 * sum_k{w_k*l_i(x_k)*[dl_j/dxi(x_k)*dxi/dx(x_k)
 *						+ dl_j/deta(x_k)*deta/dx(x_k)]*det[J](x_k)}, and
 $ sum_k{w_k*l_i(x_k)*[dl_j/dxi(x_k)*dxi/dy(x_k)
 *						+ dl_j/deta(x_k)*deta/dy(x_k)]*det[J](x_k)}, where the
 * sum is done over the Gauss points (GP). Since the components w_k*l_i(x_k) are
 * already computed in the mesh, it suffices to compute the components of the inverse
 * change of variable between the reference and physical frame (the dX/dx), to
 * calculate the sum, and store the result in a sparse matrix.
 * \param mesh The structure that contains the mesh.
 * \param Sx The Eigen::SparseMatrix in which the matrix [Sx] components will be
 * stored.
 * \param Sy The Eigen::SparseMatrix in which the matrix [Sx] components will be
 * stored.
 */
void buildS(const Mesh& mesh, Eigen::SparseMatrix<double>& Sx,
	Eigen::SparseMatrix<double>& Sy);

#endif /* buildS_hpp */
