#ifndef buildM_hpp
#define buildM_hpp

#include <Eigen/Sparse>
#include "readMesh.hpp"

void buildM(const MeshParams& meshParams, Eigen::SparseMatrix<double>& M);

#endif /* buildM_hpp */
