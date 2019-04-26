#ifndef buildMatrix_hpp
#define buildMatrix_hpp

#include "../mesh/Mesh.hpp"
#include "matrix.hpp"

/**
 * \brief Builds the matrices required for the DG-FEM
 * \param mesh The structure that contains the mesh.
 * \param matrix The structure that will contain the matrices.
 */
void buildMatrix(const Mesh& mesh, Matrix& matrix);

#endif /* buildMatrix_hpp */
