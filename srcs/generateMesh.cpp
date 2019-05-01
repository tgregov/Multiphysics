#include <gmsh.h>
#include "generateMesh.hpp"

void generateMesh(const std::string& model, const std::string& mesh, const unsigned int order){
gmsh::initialize();
gmsh::option::setNumber("General.Terminal", 1);
gmsh::open(model);
gmsh::model::mesh::generate(2);
gmsh::model::mesh::setOrder(order);
gmsh::write(mesh);
gmsh::finalize();

}