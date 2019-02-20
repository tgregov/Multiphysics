//
//  readMesh.hpp
//  

#ifndef readMesh_hpp
#define readMesh_hpp

#include <vector>

struct Mesh
{
    // structure that contains all the information of the mesh that this 
    // necessary for the DG method
    // * elementTags = [e1, e2, e3, ...]: concatenation of element tags
    //      [maybe not]
    // * nodeElements = [[n1a, n1b, n1c], [n2a, ...], ...]: node tags 
    //      corresponding to each element
    //      [end maybe not]
    // * nodeTags = [n1a, n1b, n2a, n2b, ...]: concatenation of the edges nodes,
    // 		where each node appears only once ("reduced" list of nodes)
    // * parentElement = [[e1, e2], [e1], ...]: vector such that 
    // 		parentElement[i] is associated to the i-th edge contained in 
    //		nodeTags (i.e., the node tags of the edge are (nodesTags[2*i],
    //		nodesTags[2*i+1])
    // * normalVector = [[n1x, n1y, n2x, n2y], [n2x, n2y], ...]: vector that 
    //      contains the list of normals, where the order is similar to the one
    //      of parentElement: normalVector[i] contains the normal(s) at the edge
    //      corresponding to (nodesTags[2*i], nodesTags[2*i+1]), such that there
    //      is either one normal or two normals (i.e. 2 or 4 components)
    // * jacElement = [jac1gp1, jac1gp2, jac2gp1, ...]: jacobian of each 
    //      element, evaluated at the Gauss points
    // * jacEdges = [jac1gp1, jac1gp2, jac2gp1, ...]: jacobian of each edge, 
    //      evaluated at the Gauss points
    std::vector<int> elementTags;
    //std::vector<std::vector<int>> elementTags;
    std::vector<int> nodeTags;
   	std::vector<std::vector<int>> parentElement;
    std::vector<std::vector<double>> normalVector;
    std::vector<double> jacElement;
    std::vector<double> jacEdges;
};

int posInVector(const std::vector<int>& vec, const std::vector<int>& couple);
Mesh* readMesh(int argc, char **argv);

#endif /* readMesh_hpp */
