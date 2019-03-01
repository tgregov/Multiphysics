#ifndef buildS_hpp
#define buildS_hpp

#include <Eigen/Sparse>
#include "readMesh.hpp"

void buildS(const MeshParams& meshParams, Eigen::SparseMatrix<double>& Sx, 
			Eigen::SparseMatrix<double>& Sy);

#endif /* buildS_hpp */
