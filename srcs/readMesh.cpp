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


bool readMesh(MeshParams& meshParams, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType,
                const std::string& basisFuncGradType)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName);

    // loop over the dimension i to get the maximum element dimension in the mesh
    for(unsigned short i = 0 ; i <= 3 ; ++i)
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

    // we currently only handle 2D meshes
    if(meshParams.elementDim != 2)
    {
        gmsh::logger::write("Only 2D meshes handled currently!", "error");

        gmsh::finalize();
        return false;
    }

    std::string name;
    int dim, order, numNodes;
    std::vector<double> paramCoord;
    gmsh::model::mesh::getElementProperties(meshParams.elementType, name, dim, order, numNodes, paramCoord);


    // I. BASIS FUNCTIONS & JACOBIANS

    int numComp;
    // get basis functions and their gradients
    gmsh::model::mesh::getBasisFunctions(meshParams.elementType, intScheme,
                                            basisFuncType, meshParams.intPoints,
                                            numComp, meshParams.basisFunc);

    std::cout<<std::endl<<std::endl;

    gmsh::model::mesh::getBasisFunctions(meshParams.elementType, intScheme,
                                            basisFuncGradType, meshParams.intPoints,
                                            numComp, meshParams.basisFuncGrad);

    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, meshParams.elementDim);
    int c = entities[0].second; // c is the tag of the surface

    // get the jacobians
    std::vector<double> pts;
    gmsh::model::mesh::getJacobians(meshParams.elementType, intScheme,
                                    meshParams.jacobian, meshParams.determinant, pts,
                                    c);

    // get the number of Gauss points and shape functions
    meshParams.nGP = meshParams.intPoints.size()/4;
    meshParams.nSF = meshParams.basisFunc.size()/meshParams.nGP;


    // II. NORMALS
    // [TO DO]: handle the case where there is more than one surface
    for (std::size_t i = 0; i < entities.size(); i++)
    {
        // get the tag of the current surface
        int entityTag = entities[i].second;

        // get the tags of all the elements
        std::vector<int> nodeTags;
        gmsh::model::mesh::getElementsByType(meshParams.elementType,
                                                meshParams.elementTags, nodeTags,
                                                entityTag);

        // get the barycenters of all the elements
        std::vector<double> baryCenters;
        gmsh::model::mesh::getBarycenters(meshParams.elementType, entityTag, false,
                                            true, baryCenters);
        // [Simon] When multiple entity, table of meshParams ?

        // get the nodes (in fact, tags) on the edges of the 2D elements
        // -> to do so, we search over the 2D elements of type
        // meshParams.elementType, belonging to the entity of tag entityTag
        // N.B.: there are 2 nodes per edge per element
        std::vector<int> nodes;
        gmsh::model::mesh::getElementEdgeNodes(meshParams.elementType, nodes,
                                                entityTag, true);

        int c = gmsh::model::addDiscreteEntity(1);

        // and add new 1D elements to it, for all edges
        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodes);

        std::vector<int> nodeTagsInf;
        gmsh::model::mesh::getElementsByType(eleType1D, meshParams.elementTagsInferior, nodeTagsInf, c);

        int numInferior;
        gmsh::model::mesh::getBasisFunctions(eleType1D, intScheme, basisFuncType, meshParams.intPointsInferior,
                                             numInferior, meshParams.basisFuncInferior);

        std::vector<double> jacInf, ptsInferior;
        gmsh::model::mesh::getJacobians(eleType1D, intScheme, jacInf, meshParams.determinantInferior, ptsInferior);

        meshParams.nGPInf = meshParams.intPointsInferior.size()/4;
        meshParams.nSFInf = meshParams.basisFuncInferior.size()/meshParams.nGPInf;
        meshParams.nEInf = meshParams.elementTagsInferior.size(); //Problem 4 elements 1D too much :/

        // save the number of elements
        meshParams.nE = meshParams.elementTags.size();

        // compute the number of edges per element
        double nEdgePerEl = nodes.size()/(2*meshParams.nE);

        // compute a vector of nodes per edge per element and normal per edge
        // per element => to be generalized to 1D and 3D

        meshParams.indexInFront.resize(nodeTags.size());
        for(unsigned int i = 0 ; i < nodeTags.size() ; ++i)
        {
            std::cout<<"Lol     "<<nodeTags[i]<<std::endl;
            bool founded = false;
            for(unsigned int j = 0 ; j < nodeTags.size() ; ++j)
            {
                if(j != i)
                {
                    if(nodeTags[i] == nodeTags[j])
                    {
                        meshParams.indexInFront[i]=j;
                        founded=true;
                        break;
                    }
                }
            }
            if(!founded)
                meshParams.indexInFront[i]=-1;
            std::cout<<"Smite   "<<meshParams.indexInFront[i]<<std::endl;
        }

        // loop over the elements
        for(std::size_t elm = 0; elm < meshParams.nE; elm++)
        {
            // edgeList contains, for the current element, for all its edges, the
            // associated nodes tags
            // normalList contains, for the current element, for all its edges, the
            // components of the normal
            std::vector<std::pair<int, int>> edgeList;
            std::vector<std::vector<double>> normalList;
            std::vector<double> meanNormal(2,0);
            std::vector<int> elementIndex;
            elementIndex.push_back(nEdgePerEl*elm);
            elementIndex.push_back(nEdgePerEl*elm+1);
            elementIndex.push_back(nEdgePerEl*elm+2);
            meshParams.index.push_back(elementIndex);

            // loop over the edges of the current element
            for(unsigned short k = 0 ; k < nEdgePerEl; k++)
            {

                std::vector<double> normal;
                // get the nodes of the current edge (the "2" factor comes from the
                // fact that there are each time 2 values for a single edge: the
                // associated node tags)
                edgeList.push_back(std::pair<int, int>(nodes[2*nEdgePerEl*elm + 2*k],
                                                    nodes[2*nEdgePerEl*elm + 2*k+1]));

                // get the coordinates of the first node
                std::vector<double> coord1, parametricCoord1;
                gmsh::model::mesh::getNode(edgeList[k].first, coord1,
                                            parametricCoord1);

                // get the coordinates of the second node
                std::vector<double> coord2, parametricCoord2;
                gmsh::model::mesh::getNode(edgeList[k].second, coord2,
                                            parametricCoord2);

                // compute the normal
                // if A:(x1, y1) and B:(x2, y2), then AB = (x2 - x1, y2 - y1) and a
                // normal is given by n = (y2 - y1, x1 - x2)
                // [TO DO] generalize for 1D and 3D
                double nx = coord2[1] - coord1[1];
                double ny = coord1[0] - coord2[0];
                double norm = sqrt(ny*ny + nx*nx);

                // unfortunately, nodes per edge in nodes vector are not always in
                // the same order (clockwise vs anticlockwise) => we need to check
                // the orientation
                double vx = baryCenters[3*elm] - (coord2[0] + coord1[0])/2;
                double vy = baryCenters[3*elm + 1] - (coord2[1] + coord1[1])/2;

                if(nx*vx + ny*vy > 0)
                {
                    nx = -nx;
                    ny = -ny;
                }

                // normalize the normal components
                normal.push_back(nx/norm);
                normal.push_back(ny/norm);

                // add the normal to the list of normals of the current element
                normalList.push_back(normal);
                meanNormal[0]+=normal[0]*meshParams.determinantInferior[elm*meshParams.nGPInf];
                meanNormal[1]+=normal[1]*meshParams.determinantInferior[elm*meshParams.nGPInf];  //to check
            }

            // meshParams.nodes is a vector of vector of a pair: it lists for all the
            // elements, then for all the edges the pair of node tags
            // meshParams.normals is a vector of vector of vector: it lists for all
            // the elements, then for all the edges the components of the exterior
            // normal
            meshParams.nodes.push_back(edgeList);
            meshParams.normals.push_back(normalList);
        }

        // [TO REMOVE] display the normals to check that it works
        for(std::size_t j = 0; j < meshParams.nE; j++)
        {
            std::cout<<"Element: " << meshParams.elementTags[j] << std::endl;
            for (unsigned short i = 0 ; i < meshParams.normals[j].size() ; i++ )
            {
                std::cout << "\tEdge: " << i << std::endl;
                std::cout << "\t\tNormal: (";
                for (unsigned short k = 0 ; k < meshParams.normals[j][i].size() ; k++ )
                {
                    std::cout << meshParams.normals[j][i][k];
                    if(k != meshParams.normals[j][i].size() - 1)
                    {
                        std::cout << ", ";
                    }
                }
                std::cout << ")" << std::endl;
            }
        }
    }

    gmsh::finalize();

    return true;
}
