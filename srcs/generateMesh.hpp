#ifndef generateMesh_hpp_included
#define generateMesh_hpp_included

#include <string>


/**
 * \brief generate a 2D mesh of order "order" from the .geo file "model"
 * write it in the .msh file "mesh"
 */
void generateMesh(const std::string& model, const std::string& mesh, const unsigned int order);

void generateMeshQuad(const std::string& model, const std::string& mesh, const unsigned int order);

#endif /* generateMesh_hpp */