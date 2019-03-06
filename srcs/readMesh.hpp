#ifndef readMesh_hpp
#define readMesh_hpp
#include <Eigen/Sparse>
#include <vector>
#include <string>

/**
 * \struct Mesh
 * \brief structure that contains all the information of the mesh that this
 * necessary for the DG method
 */
struct Mesh
{
    // * nodeElements = [[n1a, n1b, n1c], [n2a, ...], ...]: node tags
    //      corresponding to each element
    //      [end maybe not]
    std::vector<int> elementTags; /*!< [e1, e2, e3, ...]: concatenation of
                                    element tags [maybe not] */
    //std::vector<std::vector<int>> elementTags;
    std::vector<int> nodeTags; /*!< [n1a, n1b, n2a, n2b, ...]: concatenation of
                                  the edges nodes, where each node appears only
                                  once ("reduced" list of nodes) */

    std::vector<std::vector<int>> parentElement; /**< [[e1, e2], [e1], ...]:
                                                    vector such that
                                                    parentElement[i] is
                                                    associated to the i-th edge
                                                    contained in nodeTags (i.e.,
                                                    the node tags of the edge
                                                    are (nodesTags[2*i],
                                                    nodesTags[2*i+1]) */

    std::vector<std::vector<double>> normalVector;/**< [[n1x, n1y, n2x, n2y],
                                                    [n2x, n2y], ...]: vector
                                                    that contains the list of
                                                    normals, where the order is
                                                    similar to the one of
                                                    parentElement:
                                                    normalVector[i] contains the
                                                    normal(s) at the edge
                                                    corresponding to
                                                    (nodesTags[2*i],
                                                    nodesTags[2*i+1]), such that
                                                    there is either one normal
                                                    or two normals (i.e. 2 or 4
                                                    components) */

    std::vector<double> jacElement;/**< [jac1gp1, jac1gp2, jac2gp1, ...]:
                                      jacobian of each element, evaluated at the
                                      Gauss points*/

    std::vector<double> jacEdges;/**<  [jac1gp1, jac1gp2, jac2gp1, ...]:
                                    jacobian of each edge, evaluated at the
                                    Gauss points*/
};

/**
 * \struct MeshParams
 * \brief structure that contains all the information of the mesh that this
 * necessary for the DG method
 */
struct MeshParams
{
    std::vector<int> elementTags;  /*!< [e1, e2, e3, ...]: concatenation of
                                    element tags [maybe not] */

    std::vector<std::vector<std::pair<int, int>>> nodes; /**< vector (length same as elementTags) containing a vector
                                                              (length = number of edges per element)
                                                              of pair of nodes tags */

    std::vector<std::vector<std::vector<double>>> normals; /**< vector (length same as elementTags) containing a vector
                                                              (length = number of edges per element)
                                                              of normals (length depends on dimension */

    std::vector<std::vector<int>> index; //index[elm]=index of nodes f element elm, in the nodes vector returned somewhere in meshparams.
    std::vector<int> indexInFront; //give it an index, gives you the index of the node in front of

    std::vector<double> basisFunc;
    std::vector<double> basisFuncGrad;
    std::vector<double> determinant;
    std::vector<double> jacobian;
    std::vector<double> intPoints;

    unsigned int nGP; //maybe short
    unsigned int nSF;
    unsigned int nE;

    std::vector<int> elementTagsInferior;

    std::vector<double> basisFuncInferior;
    std::vector<double> determinantInferior;
    std::vector<double> intPointsInferior;

    unsigned int nGPInf; //maybe short
    unsigned int nSFInf;
    unsigned int nEInf;

    int elementDim;
    int elementType;
    int elementTypeInferior = -1; //type of the boundary, dummy value for dim 0

    std::vector<Eigen::SparseMatrix<double>> dM ;
    unsigned int nSigma = 3;
};

bool readMesh(MeshParams& meshParams, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType,
                const std::string& basisFuncGradType);

#endif /* readMesh_hpp */
