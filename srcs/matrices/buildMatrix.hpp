#ifndef buildMatrix_hpp_included
#define buildMatrix_hpp_included

#include "../mesh/Mesh.hpp"
#include "../utils/utils.hpp"
#include "matrix.hpp"

/**
 * \brief Builds the matrices required for the DG-FEM
 * \param mesh The structure that contains the mesh.
 * \param matrix The structure that will contain the matrices.
 */
void buildMatrix(const Mesh& mesh, Matrix& matrix,
                 const DomainDiv& domainDiv, unsigned int rank);

#endif /* buildMatrix_hpp_included */
