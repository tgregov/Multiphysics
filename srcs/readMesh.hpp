#ifndef readMesh_hpp
#define readMesh_hpp

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

int posInVector(const std::vector<int>& vec, const std::vector<int>& couple);
bool readMesh(Mesh& mesh, const std::string& fileName);

#endif /* readMesh_hpp */
