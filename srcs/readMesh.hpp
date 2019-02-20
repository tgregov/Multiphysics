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
    // * nodeTags = [n1a, n1b, n2a, n2b, ...]: concatenation of the edges nodes,
    // 		where each node appears only once ("reduced" list of nodes)
    // * parentElement = [[e1, e2], [e1], ...]: vector such that 
    // 		parentElement[i] is associated to the i-th edge contained in 
    //		nodeTags (i.e., the node tags of the edge are (nodesTags[2*i],
    //		nodesTags[2*i+1])
    std::vector<int> elementTags;
    std::vector<int> nodeTags;
   	std::vector<std::vector<int>> parentElement;
};

int posInVector(const std::vector<int>& vec, const std::vector<int>& couple);
Mesh* readMesh(int argc, char **argv);

#endif /* readMesh_hpp */
