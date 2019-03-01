#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
#include <gmsh.h>
#include <cmath>
#include "readMesh.hpp"


int posInVector(const std::vector<int>& vec, const std::pair<int, int>& couple)
{
    for(std::size_t i = 0; i < vec.size()/2; i++)
    {
        if(vec[2*i] == couple.first && vec[2*i+1] == couple.second)
            return i;
        if(vec[2*i] == couple.second && vec[2*i+1] == couple.first)
            return i;
    }

    return -1;
}


// TO DO: optimize this procedure (linked to the normal computation)
int diffVector(const std::vector<int>& vec, const std::pair<int, int>& couple)
{
    for(std::size_t i = 0; i < vec.size(); i++){
        if(vec[i] != couple.first && vec[i] != couple.second)
            return vec[i];
    }

    return -1;
}


bool readMesh(MeshParams& meshParams, const std::string& fileName, const std::string& intScheme, const std::string& basisFuncType)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName);

    for(unsigned short i = 1 ; i <= 3 ; ++i)
    {
        std::vector<int> eleTypes;
        gmsh::model::mesh::getElementTypes(eleTypes, i);
        switch(eleTypes.size())
        {
            case 0:
                break;
            case 1:
                meshParams.elementDim = i;
                meshParams.elementType = eleTypes[0]; // e.g. T3 elements
                break;
            default:
                gmsh::logger::write("Hybrid meshes not handled in this example!",
                                "error");

                gmsh::finalize();
                return false;
        }
    }

    if(meshParams.elementDim != 2)
    {
        gmsh::logger::write("Only 2D meshes handled currently!",
                                "error");

        gmsh::finalize();
        return false;
    }

    // get basis functions
    int numComp;
    gmsh::model::mesh::getBasisFunctions(meshParams.elementType, intScheme, basisFuncType,
                                         meshParams.intPoints, numComp, meshParams.basisFunc);

    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2);
    int c = entities[0].second; // c is the tag of the surface

    // Get the Jacobians information for the 2D triangular elements
    std::vector<double> jac, pts;
    gmsh::model::mesh::getJacobians(meshParams.elementType, intScheme, jac, meshParams.determinant, pts, c);

    // TO DO: what if there are more than one surface ? -> handle that case
    for (std::size_t i = 0; i < entities.size(); i++)
    {
        // get the tag of the current surface
        int entityTag = entities[i].second;

        std::vector<int> nodeTags;
        gmsh::model::mesh::getElementsByType(meshParams.elementType, meshParams.elementTags, nodeTags, entityTag);

        std::vector<double> baryCenters;
        gmsh::model::mesh::getBarycenters(meshParams.elementType, entityTag, false, true, baryCenters);//When multiple entity, table of meshParams ? ^

        // get the nodes (in fact, tags) on the edges of the 2D elements
        // -> to do so, we search over the 2D elements of type eleType2D,
        // belonging to the entity entityTag
        std::vector<int> nodes;  // 2nodes per edge per element
        gmsh::model::mesh::getElementEdgeNodes(meshParams.elementType, nodes, entityTag, true);

        meshParams.nE = meshParams.elementTags.size();

        double nEdgePerEl = nodes.size()/(2* meshParams.nE);

        // Compute a vector of nodes per edge per element and normal per edge per element
        // to be generalized to 1D and 3D
        for(std::size_t j = 0; j < meshParams.nE; j++)
        {
            std::vector<std::pair<int, int>> edgeList;
            std::vector<std::vector<double>> normalList; //List of normal for an element
            for(unsigned short k = 0 ; k < nEdgePerEl ; k++) //Normally now any element type
            {
                std::vector<double> normal; //The normal of one edge

                edgeList.push_back(std::pair<int, int>(nodes[2*nEdgePerEl*j+2*k], nodes[2*nEdgePerEl*j+2*k+1]));

                std::vector<double> coord1, parametricCoord1;
                gmsh::model::mesh::getNode(edgeList[k].first,  coord1, parametricCoord1);

                std::vector<double> coord2, parametricCoord2;
                gmsh::model::mesh::getNode(edgeList[k].second,  coord2, parametricCoord2);

                double nx = coord2[1]-coord1[1];  //Not general, 2D case, to be genralized for 1D and 3D
                double ny = coord1[0]-coord2[0];
                double norm = sqrt(ny*ny + nx*nx);

                //Unfortunately, nodes per edge in nodes vector are not always in the same order
                //clockwise vs anticlockwise
                double vx = baryCenters[3*j]-(coord2[0]+coord1[0])/2;
                double vy = baryCenters[3*j+1]-(coord2[1]+coord1[1])/2;

                if(nx*vx + ny*vy > 0)
                {
                    nx = -nx;
                    ny = -ny;
                }

                normal.push_back(nx/norm);
                normal.push_back(ny/norm);

                normalList.push_back(normal); //we add the normal the the list of element's normal
            }

            meshParams.nodes.push_back(edgeList);
            meshParams.normals.push_back(normalList); //we add the list of normal of one element to the full list of normal.
        }

        //Just read what we have found
        for(std::size_t j = 0; j < meshParams.nE; j++)
        {
            std::cout<<"Element: "<<meshParams.elementTags[j]<<std::endl;
            for (unsigned short i = 0 ; i < meshParams.normals[j].size() ; i++ )
            {
                std::cout<<"\tEdge: "<<i<<std::endl;
                std::cout<<"\t\tNormal: (";
                for (unsigned short k = 0 ; k < meshParams.normals[j][i].size() ; k++ )
                {
                    std::cout<<meshParams.normals[j][i][k]<<", ";
                }
                std::cout<<")"<<std::endl;
            }
        }
    }

    gmsh::finalize();

    meshParams.nGP = meshParams.intPoints.size()/4;
    meshParams.nSF = meshParams.basisFunc.size()/meshParams.nGP;

    return true;
}
