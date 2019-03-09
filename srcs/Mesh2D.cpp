/**
 * \file Mesh2D.cpp
 * \brief Implementation of the required function to load a Mesh2D
 * struct from file.
 *
 */

#include <iostream>
#include <gmsh.h>
#include "Mesh2D.hpp"

/**
 * \brief Loads the name order, dimension, number of nodes,
 *  basis functions, integration points for a certain element type into a map.
 * \param meshElementProp The map in which the element properties are stored.
 * \param eleTypes A vector of element types for which properties will be loaded.
 * The function does not try to load them twice if it already exists.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 * \param basisFuncGradType The type of basis function you will use ("Grad" prefix)(will be droped).
 */
static void loadElementProperties(std::map<int, ElementProperty>& meshElementProp, const std::vector<int> eleTypes,
                                  const std::string& intScheme, const std::string& basisFuncType,
                                  const std::string& basisFuncGradType)
{
    for(unsigned int i = 0 ; i < eleTypes.size() ; ++i)
    {
        if(meshElementProp.count(eleTypes[i]) == 0)
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
}

/**
 * \brief Add an edge to a certain element (filling the required fields).
 * \param element The parent element
 * \param nodesTagsEdge Node tags per edge of the element
 */
static void addEdge(Element2D& element, std::vector<int> nodesTagsEdge)
{

}

/**
 * \brief Compute the outward normal of an edge
 * \param element the parent element
 * \param nodesTagsEdge Node tags per edge of the element
 */
static void computeEdgeNormal(Element2D& element, const std::vector<int>& nodesTagsEdge)
{

}

/**
 * \brief Add an element to a certain entity (filling the required fields).
 * \param entity The parent entity.
 * \param elementTag Element tag.
 * \param eleType2D Type of the element.
 * \param eleType1D Type of the element's edges.
 * \param jacobians2D Jacobian matrix of the variable change of the element
 * evaluated at each Gauss points.
 * \param determinants2D Determinant of the variable change of the element
 * evaluated at each Gauss points.
 * \param nodesTagsPerEdge Node tags of the element, per edge.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 * \param basisFuncGradType The type of basis function you will use ("Grad" prefix)(will be droped).
 */
static void addElement(Entity2D& entity, int elementTag, int eleType2D, int eleType1D,
                       std::vector<double> jacobians2D, std::vector<double> determinants2D, std::vector<int> nodesTagsPerEdge,
                       const std::string& intScheme, const std::string& basisFuncType, const std::string& basisFuncGradType)
{
    Element2D element;
    element.elementTag = elementTag;
    element.elementType2D = eleType2D;
    element.elementType1D = eleType1D;

    element.determinant2D = determinants2D;
    element.jacobian2D = jacobians2D;

    for(unsigned int i = 0 ; i < nodesTagsPerEdge.size()/2 ; ++i)
    {
        std::vector<int> nodesTagsEdge(nodesTagsPerEdge.begin() + 2*i, nodesTagsPerEdge.begin() + 2*(i + 1));
        computeEdgeNormal(element, nodesTagsEdge);
        addEdge(element, nodesTagsEdge);
    }

    entity.elements.push_back(element);
}

/**
 * \brief Add an entity to a certain 2D mesh (filling the required fields).
 * \param mesh The parent mesh.
 * \param entityHandle Entity dimTags to add.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 * \param basisFuncGradType The type of basis function you will use ("Grad" prefix)(will be droped).
 */
static void addEntity(Mesh2D& mesh, const std::pair<int, int>& entityHandle,
                      const std::string& intScheme, const std::string& basisFuncType, const std::string& basisFuncGradType)
{
    Entity2D entity;
    entity.entityTag2D = entityHandle.second;
    int c = gmsh::model::addDiscreteEntity(1);
    entity.entityTag1D = c;

    std::vector<int> eleTypes2D;
    gmsh::model::mesh::getElementTypes(eleTypes2D, entityHandle.first, entityHandle.second);

    loadElementProperties(mesh.elementProperties2D, eleTypes2D,
                          intScheme, basisFuncType, basisFuncGradType);

    for(auto eleType2D : eleTypes2D)
    {
        std::vector<int> nodeTags, elementTags;
        gmsh::model::mesh::getElementsByType(eleType2D, elementTags, nodeTags, entityHandle.second);
        entity.elementTags2D[eleType2D] = elementTags;
        entity.nodesTags2D[eleType2D] = nodeTags;

        std::vector<int> nodesTagPerEdge;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodesTagPerEdge, entityHandle.second);
        entity.nodesTagsPerEdge2D[eleType2D] = nodesTagPerEdge;

        // Add 1D entity to store all the line associated to elements of the same order
        // Problem Q4 and T3 element have the same order ?
        int eleType1D = gmsh::model::mesh::getElementType("line", mesh.elementProperties2D[eleType2D].order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodesTagPerEdge);

        loadElementProperties(mesh.elementProperties1D, std::vector<int>(1, eleType1D), intScheme, basisFuncType, basisFuncGradType);

        std::vector<double> jacobians2D, determinants2D, dummyPoints2D;
        gmsh::model::mesh::getJacobians(eleType2D, intScheme, jacobians2D, determinants2D, dummyPoints2D, entityHandle.second);

        unsigned int nGP2D = mesh.elementProperties2D[eleType2D].intPoints.size()/4;
        unsigned int nEdgePerNode = nodesTagPerEdge.size()/nodeTags.size(); //TO check
        for(unsigned int i = 0 ; i < elementTags.size() ; ++i)
        {
            std::vector<double> jacobiansElement2D(jacobians2D.begin() + 9*nGP2D*i, jacobians2D.begin() + 9*nGP2D*(1 + i));
            std::vector<double> determinantsElement2D(determinants2D.begin() + nGP2D*i, determinants2D.begin() + nGP2D*(1 + i));
            std::vector<int> nodesTagPerEdgeElement(nodesTagPerEdge.begin() + 2*nEdgePerNode*i, nodesTagPerEdge.begin() + 2*nEdgePerNode*(i + 1));
            addElement(entity, elementTags[i], eleType2D, eleType1D, std::move(jacobiansElement2D), std::move(determinantsElement2D),
                       std::move(nodesTagPerEdgeElement), intScheme, basisFuncType, basisFuncGradType);
        }
    }

    mesh.entities.push_back(entity);
}

/**
 * \brief Checks if the mesh is 2D.
 * \return true if the mesh is 2D, false otherwise.
 */
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

// Doc in .hpp
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

    std::vector<std::pair<int, int>> entityHandles;
    gmsh::model::getEntities(entityHandles, 2);

    for(auto entityHandle : entityHandles)
    {
        addEntity(mesh, entityHandle, intScheme, basisFuncType, basisFuncGradType);
    }

    gmsh::finalize();

    return true;
}
