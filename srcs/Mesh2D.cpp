/**
 * \file Mesh2D.cpp
 * \brief Implementation of the required function to load a Mesh2D struct from file.
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
 */
static void loadElementProperties(std::map<int, ElementProperty>& meshElementProp,
                                    const std::vector<int>& eleTypes,
                                    const std::string& intScheme,
                                    const std::string& basisFuncType)
{
    for(unsigned int i = 0 ; i < eleTypes.size() ; ++i)
    {
        // if the element type has already been considered, there is no need to
        // add its properties in the map meshElementProp
        if(meshElementProp.count(eleTypes[i]) == 0)
        {
            ElementProperty elementProperty;
            gmsh::model::mesh::getElementProperties(eleTypes[i],
                                                    elementProperty.name,
                                                    elementProperty.dim,
                                                    elementProperty.order,
                                                    elementProperty.numNodes,
                                                    elementProperty.paramCoord);

            gmsh::model::mesh::getBasisFunctions(eleTypes[i], intScheme,
                                                    basisFuncType,
                                                    elementProperty.intPoints,
                                                    elementProperty.numComp,
                                                    elementProperty.basisFunc);

            std::vector<double> dummyIntPoints;
            int dummyNumComp;
            gmsh::model::mesh::getBasisFunctions(eleTypes[i], intScheme,
                                                    std::string("Grad" + basisFuncType),
                                                    dummyIntPoints,
                                                    dummyNumComp,
                                                    elementProperty.basisFuncGrad);


            // add the products prodFunc[k][i,j] = w_k*l_i(x_k)*l_j(x_k)
            // add the products pondFunc[k][i] = w_k*l_i(x_k)
            elementProperty.nGP = elementProperty.intPoints.size()/4;
            elementProperty.nSF = elementProperty.basisFunc.size()/elementProperty.nGP;
            for(unsigned int k = 0 ; k < elementProperty.nGP ; ++k)
            {
                std::vector<double> wll;
                std::vector<double> wl;

                for(unsigned int i = 0 ; i < elementProperty.nSF ; ++i)
                {

                    wl.push_back(elementProperty.intPoints[4*k + 3]
                            *elementProperty.basisFunc[elementProperty.nSF*k + i]);

                    for(unsigned int j = i ; j < elementProperty.nSF ; ++j)
                    {
                        if(k == 0)
                        {
                            elementProperty.IJ.push_back(
                                std::pair<unsigned int, unsigned int>(i, j));
                        }

                        wll.push_back(elementProperty.intPoints[4*k + 3]
                            *elementProperty.basisFunc[elementProperty.nSF*k + i]
                            *elementProperty.basisFunc[elementProperty.nSF*k + j]);
                    }
                }

                elementProperty.pondFunc.push_back(wl);
                elementProperty.prodFunc.push_back(wll);
            }


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

    Edge edge;
    edge.nodeTags.first = nodesTagsEdge[0];
    edge.nodeTags.second = nodesTagsEdge[1];
    element.edges.push_back(edge);
    // we still need to add its tag & det

}


/**
 * \brief Compute the outward normal of an edge
 * \param element the parent element
 * \param nodesTagsEdge Node tags per edge of the element
 * \param baryCenter  Barycenter of the parent element
 */
static void computeEdgeNormal(Element2D& element,
                                const std::vector<int>& nodesTagsEdge,
                                const std::vector<double>& baryCenter)
{
    std::pair<double, double> normal;

    std::vector<double> coord1, dummyParametricCoord1;
    gmsh::model::mesh::getNode(nodesTagsEdge[0], coord1,
                                            dummyParametricCoord1);

    // get the coordinates of the second node
    std::vector<double> coord2, dummyParametricCoord2;
    gmsh::model::mesh::getNode(nodesTagsEdge[1], coord2,
                                dummyParametricCoord2);

    // compute the normal
    // if A:(x1, y1) and B:(x2, y2), then AB = (x2 - x1, y2 - y1) and a
    // normal is given by n = (y2 - y1, x1 - x2)
    double nx = coord2[1] - coord1[1];
    double ny = coord1[0] - coord2[0];
    double norm = sqrt(ny*ny + nx*nx);

    // unfortunately, nodes per edge in nodes vector are not always in
    // the same order (clockwise vs anticlockwise) => we need to check
    // the orientation
    double vx = baryCenter[0] - (coord2[0] + coord1[0])/2;
    double vy = baryCenter[1] - (coord2[1] + coord1[1])/2;

    if(nx*vx + ny*vy > 0)
    {
        nx = -nx;
        ny = -ny;
    }

    // normalize the normal components
    normal.first = nx/norm;
    normal.second = ny/norm;

    element.edgesNormal.push_back(normal);
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
 */
static void addElement(Entity2D& entity, int elementTag, int eleType2D,
                        int eleType1D, std::vector<double> jacobians2D,
                        std::vector<double> determinants2D,
                        const std::vector<int>& nodesTagsPerEdge,
                        const std::vector<double>& elementBarycenter,
                        const std::string& intScheme,
                        const std::string& basisFuncType)
{
    Element2D element;
    element.elementTag = elementTag;
    element.elementType2D = eleType2D;
    element.elementType1D = eleType1D;

    element.determinant2D = std::move(determinants2D);
    element.jacobian2D = std::move(jacobians2D);

    for(unsigned int i = 0 ; i < nodesTagsPerEdge.size()/2 ; ++i)
    {
        std::vector<int> nodesTagsEdge(nodesTagsPerEdge.begin() + 2*i,
                                        nodesTagsPerEdge.begin() + 2*(i + 1));

        computeEdgeNormal(element, nodesTagsEdge, elementBarycenter);
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
 */
static void addEntity(Mesh2D& mesh, const std::pair<int, int>& entityHandle,
                      const std::string& intScheme, const std::string& basisFuncType)
{
    // add the current 2D entity
    Entity2D entity;
    entity.entityTag2D = entityHandle.second;
    int c = gmsh::model::addDiscreteEntity(1);
    entity.entityTag1D = c;

    // get the element types in the current 2D entity
    std::vector<int> eleTypes2D;
    gmsh::model::mesh::getElementTypes(eleTypes2D, entityHandle.first,
                                        entityHandle.second);
    loadElementProperties(mesh.elementProperties2D, eleTypes2D,
                          intScheme, basisFuncType);

    // loop over the element types in the current 2D entity
    for(auto eleType2D : eleTypes2D)
    {
        // get the elements of type eleType2D
        std::vector<int> nodeTags, elementTags;
        gmsh::model::mesh::getElementsByType(eleType2D, elementTags, nodeTags,
                                                entityHandle.second);
        entity.elementTags2D[eleType2D] = elementTags;
        entity.nodesTags2D[eleType2D] = nodeTags;

        std::vector<int> nodesTagPerEdge;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodesTagPerEdge,
                                                entityHandle.second);
        entity.nodesTagsPerEdge2D[eleType2D] = nodesTagPerEdge;

        std::vector<double> baryCenters;
        gmsh::model::mesh::getBarycenters(eleType2D, entityHandle.second, false, true, baryCenters);

        // add 1D entity to store all the lines associated to elements of the same
        // order
        // TO DO: Problem Q4 and T3 element have the same order ?
        int eleType1D = gmsh::model::mesh::getElementType("line",
                                        mesh.elementProperties2D[eleType2D].order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodesTagPerEdge);

        loadElementProperties(mesh.elementProperties1D,
                                std::vector<int>(1, eleType1D), intScheme,
                                basisFuncType);

        std::vector<double> jacobians2D, determinants2D, dummyPoints2D;
        gmsh::model::mesh::getJacobians(eleType2D, intScheme, jacobians2D,
                                        determinants2D, dummyPoints2D,
                                        entityHandle.second);

        unsigned int nGP2D = mesh.elementProperties2D[eleType2D].intPoints.size()/4;
        // TO DO: check
        // unsigned int nEdgePerNode = nodesTagPerEdge.size()/nodeTags.size();
        // std::cout << " ====> " << nEdgePerNode << std::endl;
        unsigned int numNodes = mesh.elementProperties2D[eleType2D].numNodes;
        // std:: cout << numNodes << std::endl;

        for(unsigned int i = 0 ; i < elementTags.size() ; ++i)
        {
            std::vector<double> jacobiansElement2D(
                                            jacobians2D.begin() + 9*nGP2D*i,
                                            jacobians2D.begin() + 9*nGP2D*(1 + i));
            std::vector<double> determinantsElement2D(
                                            determinants2D.begin() + nGP2D*i,
                                            determinants2D.begin() + nGP2D*(1 + i));
            std::vector<int> nodesTagPerEdgeElement(
                                nodesTagPerEdge.begin() + 2*numNodes*i,
                                nodesTagPerEdge.begin() + 2*numNodes*(i + 1));

            std::vector<double> elementBarycenter(baryCenters.begin() + 3*i, baryCenters.begin() + 3*(i + 1));

            addElement(entity, elementTags[i], eleType2D, eleType1D,
                        std::move(jacobiansElement2D),
                        std::move(determinantsElement2D),
                        nodesTagPerEdgeElement,
                        elementBarycenter, intScheme, basisFuncType);
        }
    }

    // add the entity to the mesh.entities field
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
                std::cerr   << "Hybrid meshes not handled in this example!"
                            << std::endl;
                return false;
        }
    }

    if(elementDim != 2)
    {
        std::cerr   << "The mesh does not seem to be 2D "
                    << "(contains elements of higher dimension)" << std::endl;
        return false;
    }

    return true;
}


// documentation in .hpp file
unsigned long getNumNodes(const Mesh2D& mesh2D)
{
    unsigned long numNodes = 0;

    // loop over the entities
    for(unsigned int ent = 0 ; ent < mesh2D.entities.size() ; ++ent)
    {
        Entity2D entity = mesh2D.entities[ent];

        // loop over the elements
        for(unsigned int elm = 0 ; elm < entity.elements.size() ; ++elm)
        {
            // the number of nodes for an element equals its number of edges
            numNodes += entity.elements[elm].edges.size();
        }
    }

    return numNodes;
}



// documentation in .hpp file
bool readMesh2D(Mesh2D& mesh, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName);

    // check that the mesh is 2D
    if(!IsMesh2D())
    {
        gmsh::finalize();
        return false;
    }

    // collect the information contained in the gmsh file
    std::vector<std::pair<int, int>> entityHandles;
    gmsh::model::getEntities(entityHandles, 2);

    for(auto entityHandle : entityHandles)
    {
        addEntity(mesh, entityHandle, intScheme, basisFuncType);
    }

    gmsh::finalize();

    return true;
}
