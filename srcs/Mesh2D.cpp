#include <iostream>
#include <gmsh.h>
#include "Mesh2D.hpp"

static bool IsMesh2D()
{
    int elementDim = -1;

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
                elementDim = i;
                break;
            default:
                elementDim = i;
                std::cerr<<"Hybrid meshes not handled in this example!"<<std::endl;
                return false;
        }
    }

    if(elementDim != 2)
    {
        std::cerr<<"The mesh does not seem to be 2D (contains elements of higher dimension)"<<std::endl;
        return false;
    }

    return true;
}

static void loadElementProperties(std::map<int, ElementProperty> meshElementProp, const std::vector<int> eleTypes,
                                  const std::string& intScheme, const std::string& basisFuncType,
                                  const std::string& basisFuncGradType)
{
    for(unsigned int i = 0 ; i < eleTypes.size() ; ++i)
    {
        ElementProperty elementProperty;
        gmsh::model::mesh::getElementProperties(eleTypes[i], elementProperty.name, elementProperty.dim,
                                                elementProperty.order, elementProperty.numNodes,
                                                elementProperty.paramCoord);

        gmsh::model::mesh::getBasisFunctions(eleTypes[i], intScheme,
                                            basisFuncType, elementProperty.intPoints,
                                            elementProperty.numComp, elementProperty.basisFunc);

        std::vector<double> dummyIntPoints;
        int dummyNumComp;
        gmsh::model::mesh::getBasisFunctions(eleTypes[i], intScheme,
                                                basisFuncGradType, dummyIntPoints,
                                                dummyNumComp, elementProperty.basisFuncGrad);

        meshElementProp[eleTypes[i]] = elementProperty;
    }
}

static void addEntity(Mesh2D mesh, const std::pair<int, int>& entityHandle, const std::vector<int> eleTypes2D,
                      const std::vector<int> eleTypes1D, const std::string& intScheme,
                      const std::string& basisFuncType, const std::string& basisFuncGradType)
{
    Entity2D entity;
    gmsh::model::mesh::getElementTypes(eleTypes2D, 2, entityHandle.second);
    for(unsigned int i = 0 ; i < eleTypes2D.size() ; ++i)
    {
        std::vector<double> jacobians, determinants, dummyPoints;
        gmsh::model::mesh::getJacobians(eleTypes2D[i], intScheme, jacobians, determinants, dummyPoints, entityHandle.second);
        entity.jacobian2D[eleTypes2D[i]] = jacobians;
        entity.determinant2D[eleTypes2D[i]] = determinants;
    }

    for(unsigned int i = 0 ; i < eleTypes1D.size() ; ++i)
    {
        std::vector<double> jacobians, determinants, dummyPoints;
        gmsh::model::mesh::getJacobians(eleTypes1D[i], intScheme, jacobians, determinants, dummyPoints, entityHandle.second);
        entity.jacobian1D[eleTypes1D[i]] = jacobians;
        entity.determinant1D[eleTypes1D[i]] = determinants;
    }

    std::vector<int> nodeTags;
    gmsh::model::mesh::getElementsByType(meshParams.elementType,
                                            meshParams.elementTags, nodeTags,
                                            entityTag);

    std::vector<double> baryCenters;
    gmsh::model::mesh::getBarycenters(meshParams.elementType, entityTag, false,
                                        true, baryCenters);

    std::vector<int> nodes;
    gmsh::model::mesh::getElementEdgeNodes(meshParams.elementType, nodes,
                                            entityTag, true);

    int c = gmsh::model::addDiscreteEntity(1);

    // and add new 1D elements to it, for all edges // this works only for edge of the same order !!!
    int eleType1D = gmsh::model::mesh::getElementType("line", order);
    gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodes);

    mesh.entities.push_back(entity);
}

bool readMesh2D(Mesh2D& mesh, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType,
                const std::string& basisFuncGradType)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName);

    if(!IsMesh2D())
    {
        gmsh::finalize();
        return false;
    }

    std::vector<int> eleTypes2D, eleTypes1D;
    gmsh::model::mesh::getElementTypes(eleTypes2D, 2);
    gmsh::model::mesh::getElementTypes(eleTypes1D, 1);

    loadElementProperties(mesh.elementProperties1D, eleTypes1D,
                          intScheme, basisFuncType, basisFuncGradType);

    loadElementProperties(mesh.elementProperties2D, eleTypes2D,
                          intScheme, basisFuncType, basisFuncGradType);

    std::vector<std::pair<int, int>> entityHandles;
    gmsh::model::getEntities(entityHandles, 2);

    for(auto entityHandle : entityHandles)
    {
        addEntity(mesh, entityHandle, eleTypes2D, eleTypes1D, intScheme,
                  basisFuncType, basisFuncGradType);
    }

    gmsh::finalize();

    return true;
}
