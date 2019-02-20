#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI
#include <gmsh.h>
#include "readMesh.hpp"


int posInVector(const std::vector<int>& vec, const std::vector<int>& couple)
{
    for(std::size_t i = 0; i < vec.size()/2; i++)
    {
        if(vec[2*i] == couple[0] && vec[2*i+1] == couple[1])
            return i;
        if(vec[2*i] == couple[1] && vec[2*i+1] == couple[0])
            return i;
    }
    
    return -1;
}


Mesh* readMesh(int argc, char **argv)
{
    
    // check that a .msh file was introduced
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return NULL;
    } 
        
    // initialize gmsh, open a terminal and open the file
    Mesh meshInit;
    Mesh* mesh = &meshInit;

    gmsh::initialize(argc, argv);
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(argv[1]);
    
    // get the properties of 2D elements
    std::vector<int> eleTypes;
    gmsh::model::mesh::getElementTypes(eleTypes, 2);
    if (eleTypes.size() != 1)
    {
        // TO DO: handle hybrid meshes
        gmsh::logger::write("Hybrid meshes not handled in this example!",
                            "error");
        return NULL;
    }
    int eleType2D = eleTypes[0];
    std::string name;
    int dim, order, numNodes;
    std::vector<double> paramCoord;
    gmsh::model::mesh::getElementProperties(eleType2D, name, dim, order,
                                            numNodes, paramCoord);
    
    // create a pair vector that will contain the (dim, tag) of the geometrical
    // entities, and get the entities over dim = 2, that is, the surfaces
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2);
    
    // loop over the surfaces
    // TO DO: what if there are more than one surface ? -> handle that case
    for (std::size_t i = 0; i < entities.size(); i++)
    {
        // get the tag of the current surface
        int entityTag = entities[i].second;
        gmsh::logger::write("[Surface " + std::to_string(entityTag)  + "]");
        
        // get the elements (their tag & node tag) of the current surface
        // -> to do so, we search over the 2D elements of type eleType2D,
        // belonging to the entity entityTag
        std::vector<int> elementTags, nodeTags;
        gmsh::model::mesh::getElementsByType(eleType2D, elementTags, nodeTags,
                                             entityTag);
        gmsh::logger::write("------------------------------------------------");
        gmsh::logger::write("-> there are " + std::to_string(elementTags.size())
                            + " elements:");
        for(std::size_t j = 0; j < elementTags.size(); j++)
        {
            gmsh::logger::write("   - element [" + std::to_string(j) +
                                "] has tag: " + std::to_string(elementTags[j]));
        }
        
        // get the nodes (in fact, tags) on the edges of the 2D elements
        // -> to do so, we search over the 2D elements of type eleType2D,
        // belonging to the entity entityTag
        std::vector<int> nodes;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodes, entityTag,
                                               true);
        gmsh::logger::write("------------------------------------------------");
        gmsh::logger::write("-> there are " + std::to_string(nodes.size())
                            + " (non-unique) nodes");
        gmsh::logger::write("   check: " + std::to_string(nodes.size())
                            + " = 6*"
                            + std::to_string(elementTags.size())
                            + " => ok");
        
        // reduce the number of edges in order not to have any duplicate
        // edges, and store the parent element(s) of each edges
        // TO DO:
        //  - optimize the code
        //  - maybe store the results in a structure
        std::vector<std::vector<int>> parentElement;
        std::vector<int> nodesReduced;
        int p;
        for(std::size_t j = 0; j < elementTags.size(); j++)
        {
            for(int k = 0; k < 3; k++)
            {
                std::vector<int> currentEdge;
                currentEdge.push_back(nodes[6*j + 2*k]);
                currentEdge.push_back(nodes[6*j + 2*k + 1]);
                
                p = posInVector(nodesReduced, currentEdge);
                if(p == -1){
                    nodesReduced.push_back(currentEdge[0]);
                    nodesReduced.push_back(currentEdge[1]);
                    
                    std::vector<int> temp;
                    temp.push_back(elementTags[j]);
                    parentElement.push_back(temp);
                } else{
                    parentElement[p].push_back(elementTags[j]);
                }
            }
        }
        
        gmsh::logger::write("------------------------------------------------");
        gmsh::logger::write("-> there are "
                            + std::to_string(nodesReduced.size())
                            + " reduced nodes => "
                            + std::to_string(nodesReduced.size()/2)
                            + " unique edges");
        for(std::size_t j = 0; j < parentElement.size(); j++)
        {
            gmsh::logger::write("   - edge {" + std::to_string(j) + "}: ("
                                + std::to_string(nodesReduced[2*j])
                                + ")-("
                                + std::to_string(nodesReduced[2*j+1])
                                + ") has "
                                + std::to_string(parentElement[j].size())
                                + " parent element(s):");
            for(std::size_t k = 0; k < parentElement[j].size(); k++){
                gmsh::logger::write("      - element ["
                                    + std::to_string(parentElement[j][k])
                                    + "]");
            }
        }
        gmsh::logger::write("------------------------------------------------");
        
        // create a new discrete entity of dimension 1
        int disEntity = gmsh::model::addDiscreteEntity(1);
        
        // and add new 1D elements to it, for all edges
        // N.B.: setElementsByType sets the elements of type eleType1D (here
        // lines) in the entity of dimension 1 (which is the dimension) and tag
        // disEntity by giving them the tags contained in nodesReduced
        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, disEntity, eleType1D, {},
                                             nodesReduced);


        // Save the data in a mesh structure
        // TO DO: handle the case of multiple surfaces
        mesh->elementTags = elementTags;
        mesh->nodeTags =  nodesReduced;
        mesh->parentElement = parentElement;
    }
    
    // write the new .msh file
    // gmsh::write("new_mesh.msh");
    
    // basis functions for a 1D element
    gmsh::model::mesh::getElementTypes(eleTypes, 1);
    int eleType1D = eleTypes[0];
    std::vector<double> intpts, bf;
    int numComp;
    gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss3", "IsoParametric",
                                         intpts, numComp, bf);
    
    // iterate over the new 1D elements and get integration information
    // TO DO: dynamically take into account the fact that there might be a
    // vector of entities disEntity (here I just took the last entity of dim =
    // 1) (might be useless if the domain in 2D is simple, e.g. connex)
    gmsh::model::getEntities(entities, 1);
    int c = entities[entities.size()-1].second;
    std::vector<int> elementTags, nodeTags;
    gmsh::model::mesh::getElementsByType(eleType1D, elementTags, nodeTags, c);
    std::vector<double> jac, det, pts;
    gmsh::model::mesh::getJacobians(eleType1D, "Gauss3", jac, det, pts, c);
    
    // TO DO: see what it means
    //gmsh::fltk::run();
    
    // finalize gmsh
    gmsh::finalize();


    return mesh;
}
