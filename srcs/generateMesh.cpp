#include <gmsh.h>
#include "generateMesh.hpp"
#include<iostream>

void generateMesh(const std::string& model, const std::string& mesh, const unsigned int order){
gmsh::initialize();
gmsh::option::setNumber("General.Terminal", 1);
gmsh::open(model);
gmsh::model::mesh::generate(2);
gmsh::model::mesh::setOrder(order);
gmsh::write(mesh);
gmsh::finalize();

}

void generateMeshQuad(const std::string& model, const std::string& mesh, const unsigned int order){
gmsh::initialize();
gmsh::option::setNumber("General.Terminal", 1);
gmsh::open(model);
gmsh::vectorpair dimTags;
gmsh::model::getEntities(dimTags,2);
gmsh::model::mesh::setRecombine(dimTags[0].first,dimTags[0].second);
gmsh::model::mesh::generate(2);
gmsh::model::mesh::refine();
gmsh::model::mesh::setOrder(order);
gmsh::write(mesh);
gmsh::finalize();

}