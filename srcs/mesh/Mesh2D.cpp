/**
 * \file Mesh2D.cpp
 * \brief Implementation of the required function to load a Mesh2D struct from file.
 */

#include <algorithm>
#include <iostream>
#include <gmsh.h>
#include "Mesh2D.hpp"
#include "../utils.hpp"

static void loadNodeData(Mesh2D& mesh2D)
{

    unsigned int numNodes = 0;
    std::vector<int> nodeTags;
    std::vector<std::vector<double>> coord;

    // loop over the entities
    for(size_t ent = 0 ; ent < mesh2D.entities.size() ; ++ent)
    {
        Entity2D entity = mesh2D.entities[ent];

        // loop over the elements
        for(size_t elm = 0 ; elm < entity.elements.size() ; ++elm)
        {
            Element2D element = entity.elements[elm];

            for(size_t s = 0 ; s < element.edges.size() ; ++s)
            {
                Edge edge = element.edges[s];

                for(size_t n = 0 ; n < edge.nodeCoordinate.size() ; ++n)
                {
                    std::vector<double> temp;
                    temp.push_back(edge.nodeCoordinate[n].first);
                    temp.push_back(edge.nodeCoordinate[n].second);                    
                    coord.push_back(temp);

                    nodeTags.push_back(edge.nodeTags[n]);
                    numNodes++;
                }
            }
        }
    }

    mesh2D.nodeData.numNodes = numNodes;
    mesh2D.nodeData.nodeTags = nodeTags;
    mesh2D.nodeData.coord = coord;
}



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

            std::vector<std::vector<double>> lalb(elementProperty.nSF);
            for(size_t l = 0 ; l < lalb.size() ; ++l)
            {
                lalb[l].resize(elementProperty.nSF);
            }

            for(unsigned int l = 0 ; l < elementProperty.nSF*(elementProperty.nSF+1)/2 ; ++l)
            {
                double  sum = 0.0;
                for (unsigned int k = 0 ; k < elementProperty.nGP ; ++k)
                {
                    sum += elementProperty.prodFunc[k][l];
                }

                lalb[elementProperty.IJ[l].first][elementProperty.IJ[l].second] = sum;
                if(elementProperty.IJ[l].first != elementProperty.IJ[l].second)
                {
                    lalb[elementProperty.IJ[l].second][elementProperty.IJ[l].first] = sum;
                }
            }

            elementProperty.lalb = std::move(lalb);

            meshElementProp[eleTypes[i]] = std::move(elementProperty);
        }
    }
}


/**
 * \brief Add an edge to a certain element (filling the required fields).
 * \param element The parent element.
 * \param nodesTagsEdge Node tags per edge of the element.
 * \param determinant1D Determinant associated with the edge of the element.
 * \param nNodesElement Number of nodes of the parent element.
 */
static void addEdge(Element2D& element, std::vector<int> nodesTagsEdge,
                    std::vector<double> determinant1D,
                    const std::vector<int>& nodesTagsElement)
{

    Edge edge;
    //Should be before movement !
    for(unsigned int i = 0 ; i < nodesTagsEdge.size() ; ++i)
    {
        for(unsigned int j = 0 ; j < nodesTagsElement.size() ; ++j)
        {
            if(nodesTagsEdge[i] == nodesTagsElement[j]) // We should always find it!
                edge.offsetInElm.push_back(j);

        }
    }
    edge.nodeTags = std::move(nodesTagsEdge);
    edge.determinant1D = std::move(determinant1D);

    element.edges.push_back(edge);
    // we still need to add its tag
}


/**
 * \brief Compute the outward normal of an edge.
 * \param edge The edge of which the normal is computed.
 * \param baryCenter  Barycenter of the parent element.
 */
static void computeEdgeNormalCoord(Edge& edge,
                              const std::vector<double>& baryCenter)
{
    std::pair<double, double> normal;

    std::vector<double> coord1, coord2;
    for(unsigned int i = 0 ; i < edge.nodeTags.size() ; ++i)
    {
            std::vector<double> coord, dummyParametricCoord;
            gmsh::model::mesh::getNode(edge.nodeTags[i], coord,
                                        dummyParametricCoord);
            edge.nodeCoordinate.push_back(std::pair<double, double>(coord[0], coord[1]));
            //Normally first and second node tags are vertex
            if(i == 0)
                coord1 = std::move(coord);
            else if(i == 1)
                coord2 = std::move(coord);
    }

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

    edge.normal = std::move(normal);
}

/**
 * \brief Find the second position of the edge in the entity.
 * \param entity The parent entity.
 * \param currentEdge The edge of which you want to find the neighbor.
 * \param edgePos Index of the edge in its element.
 */
static void findInFrontEdge(Entity2D& entity, Edge& currentEdge, unsigned int edgePos)
{
    unsigned int elVecSize = entity.elements.size();
    unsigned int nEdgePerEl = entity.elements[0].edges.size();
    for(unsigned int elm = 0 ; elm < elVecSize ; ++elm)
    {
        for(unsigned int k = 0 ; k < nEdgePerEl ; ++k)
        {
            std::vector<int> nodesTagsInfront = entity.elements[elm].edges[k].nodeTags;
            std::vector<unsigned int> permutation1, permutation2;

            if(isPermutation(currentEdge.nodeTags, nodesTagsInfront, permutation1, permutation2))
            {
                currentEdge.edgeInFront =  std::pair<unsigned int, unsigned int>(elm, k);
                entity.elements[elm].edges[k].edgeInFront = std::pair<unsigned int, unsigned int>(elVecSize, edgePos);
                currentEdge.nodeIndexEdgeInFront = std::move(permutation1);
                entity.elements[elm].edges[k].nodeIndexEdgeInFront = std::move(permutation2);
                return;
            }
        }
    }
}

/**
 * \brief Check if an edge lies on a boundary (and which one)
 * \param nodesTagBoundaries Map which stores the tags of the nodes belonging to a certain boundary.
 * \param edge The edge which we check if it is a boundary.
 */
static bool IsBounbdary(const std::map<std::string, std::vector<int>>& nodesTagBoundaries, Edge& edge)
{
    for(std::pair<std::string, std::vector<int>> nodeTagBoundary : nodesTagBoundaries)
    {
        for(unsigned int i = 0 ; i < nodesTagBoundaries.size()/2 ; ++i)
        {
            if(std::count(nodeTagBoundary.second.begin(), nodeTagBoundary.second.end(), edge.nodeTags[0])
               && std::count(nodeTagBoundary.second.begin(), nodeTagBoundary.second.end(), edge.nodeTags[1]))
            {
                edge.bcName=nodeTagBoundary.first;

                return true;
            }
        }
    }
    return false;
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
 * \param determinants1D Determinant of the variable change of the element's edges
 * evaluated at each Gauss points.
 * \param nGP1D Number of Gauss point for integration.
 * \param offsetInU Offset of the element in the u unknown vector.
 * \param nodesTagsPerEdge Node tags of the element, per edge.
 * \param elementBarycenter Element barycenter coordinate.
 * \param nodesTagBoundaries Map which stores the tags of the nodes belonging to a certain boundary.
 * \param element2DProperty Structure containing informations about a certain element type.
 */
static void addElement(Entity2D& entity, int elementTag, int eleType2D,
                        int eleType1D, std::vector<double> jacobians2D,
                        std::vector<double> determinants2D,
                        std::vector<double> determinants1D,
                        unsigned int nGP1D, unsigned int offsetInU,
                        std::vector<int> nodesTagsPerEdge,
                        std::vector<int> nodesTags,
                        const std::vector<double>& elementBarycenter,
                        const std::map<std::string, std::vector<int>>& nodesTagBoundaries,
                        const ElementProperty& element2DProperty,
                        const ElementProperty& element1DProperty)
{
    Element2D element;
    element.elementTag = elementTag;
    element.elementType2D = eleType2D;
    element.elementType1D = eleType1D;

    element.offsetInU = offsetInU;

    element.determinant2D = std::move(determinants2D);
    element.jacobian2D = std::move(jacobians2D);
    element.nodeTags = std::move(nodesTags);

    unsigned int nEdge = element2DProperty.numNodes/element2DProperty.order;
    unsigned int nNodePerEdge = element2DProperty.order + 1;

    for(unsigned int i = 0 ; i < nEdge ; ++i)
    {
        std::vector<int> nodesTagsEdge(nodesTagsPerEdge.begin() + nNodePerEdge*i,
                                        nodesTagsPerEdge.begin() + nNodePerEdge*(i + 1));

        std::vector<double> determinantsEdge1D(determinants1D.begin() + nGP1D*i, determinants1D.begin() + nGP1D*(i + 1));

        addEdge(element, std::move(nodesTagsEdge), std::move(determinantsEdge1D), element.nodeTags);
        if(entity.elements.size() != 0)
        {
            if(!IsBounbdary(nodesTagBoundaries, element.edges[i]))
                findInFrontEdge(entity, element.edges[i], i);

        }
        else
        {
            IsBounbdary(nodesTagBoundaries, element.edges[i]);
        }

        computeEdgeNormalCoord(element.edges[i], elementBarycenter);

        Eigen::SparseMatrix<double> dMs(element2DProperty.nSF, element2DProperty.nSF);
        std::vector<Eigen::Triplet<double>> indices;

        for(size_t nA = 0 ; nA <  element.edges[i].nodeTags.size() ; ++nA)
        {
            for(size_t nB = 0 ; nB <  element.edges[i].nodeTags.size() ; ++nB)
            {
                indices.push_back(Eigen::Triplet<double>
                    (element.edges[i].offsetInElm[nA], element.edges[i].offsetInElm[nB], element1DProperty.lalb[nA][nB]));
            }
        }

        dMs.setFromTriplets(indices.begin(), indices.end());
        element.dM.push_back(dMs);
    }

    entity.elements.push_back(element);
}


/**
 * \brief Add an entity to a certain 2D mesh (filling the required fields).
 * \param mesh The parent mesh.
 * \param entityHandle Entity dimTags to add.
 * \param currentOffset Offset in u unknown vector.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 */
static bool addEntity(Mesh2D& mesh, const std::pair<int, int>& entityHandle, unsigned int& currentOffset,
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

        std::vector<int> nodesTagPerEdge;
        gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodesTagPerEdge,
                                                entityHandle.second);

        unsigned int numNodes = mesh.elementProperties2D[eleType2D].numNodes;
        unsigned int order = mesh.elementProperties2D[eleType2D].order;

        std::vector<double> baryCenters;
        gmsh::model::mesh::getBarycenters(eleType2D, entityHandle.second, false, true, baryCenters);

        // add 1D entity to store all the lines associated to elements of the same
        // order
        // TO DO: Problem Q4 and T3 element have the same order ?
        int eleType1D = gmsh::model::mesh::getElementType("line", order);
        gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodesTagPerEdge);

        loadElementProperties(mesh.elementProperties1D,
                                std::vector<int>(1, eleType1D), intScheme,
                                basisFuncType);

        std::vector<double> jacobians2D, determinants2D, dummyPoints2D;
        gmsh::model::mesh::getJacobians(eleType2D, intScheme, jacobians2D,
                                        determinants2D, dummyPoints2D,
                                        entityHandle.second);

        //Attention: assume one eleType2D per entity !! (Same eleType1D for T3 and Q4 -> bug)
        std::vector<double> dummyJacobians1D, determinants1D, dummyPoints1D;
        gmsh::model::mesh::getJacobians(eleType1D, intScheme, dummyJacobians1D,
                                        determinants1D, dummyPoints1D,
                                        c);

        unsigned int nElements = nodeTags.size()/mesh.elementProperties2D[eleType2D].numNodes;
        unsigned int nGP2D = mesh.elementProperties2D[eleType2D].nGP;
        unsigned int nGP1D = mesh.elementProperties1D[eleType1D].nGP;
        unsigned int nEdgePerElement = determinants1D.size()/(nGP1D*nElements);
        unsigned int nNodesPerEdge = mesh.elementProperties1D[eleType1D].numNodes;

        //Check for bubble nodes in eleType2D
        if(numNodes != order * nEdgePerElement)
        {
            std::cerr << "Currently unsupported element type " << mesh.elementProperties2D[eleType2D].name
                      << " inside mesh!" << std::endl;

            return false;
        }

        unsigned ratio, currentDecade = 0;
        for(unsigned int i = 0 ; i < elementTags.size() ; ++i)
        {

            // display progress
            ratio = int(100*double(i)/double(elementTags.size()));
            if(ratio >= currentDecade)
            {
                std::cout   << "\r" << "Entity [" << entity.entityTag2D << "]: "
                            << ratio << "% of the elements computed"
                            << std::flush;
                currentDecade = ratio + 1;
            }

            std::vector<double> jacobiansElement2D(
                                            jacobians2D.begin() + 9*nGP2D*i,
                                            jacobians2D.begin() + 9*nGP2D*(1 + i));
            std::vector<double> determinantsElement2D(
                                            determinants2D.begin() + nGP2D*i,
                                            determinants2D.begin() + nGP2D*(1 + i));
            std::vector<double> determinantElement1D(
                                            determinants1D.begin() + nEdgePerElement*nGP1D*i,
                                            determinants1D.begin() + nEdgePerElement*nGP1D*(1 + i));

            std::vector<int> nodeTagsElement(nodeTags.begin() + numNodes*i,
                                             nodeTags.begin() + numNodes*(1 + i));

            //[TO DO]: generalize for non-linear element
            std::vector<int> nodesTagPerEdgeElement(
                                nodesTagPerEdge.begin() + nNodesPerEdge*nEdgePerElement*i,
                                nodesTagPerEdge.begin() + nNodesPerEdge*nEdgePerElement*(i + 1));

            std::vector<double> elementBarycenter(baryCenters.begin() + 3*i, baryCenters.begin() + 3*(i + 1));

            unsigned int elementOffset = numNodes; //To check

            addElement(entity, elementTags[i], eleType2D, eleType1D,
                        std::move(jacobiansElement2D),
                        std::move(determinantsElement2D),
                        std::move(determinantElement1D),
                        nGP1D, currentOffset,
                        std::move(nodesTagPerEdgeElement),
                        std::move(nodeTagsElement),
                        elementBarycenter,
                        mesh.nodesTagBoundary,
                        mesh.elementProperties2D[eleType2D],
                        mesh.elementProperties1D[eleType1D]);

            currentOffset += elementOffset;
        }

         std::cout  << "\r" << "Entity [" << entity.entityTag2D << "]: "
                    << "100% of the elements computed" << std::flush << std::endl;
    }

    // add the entity to the mesh.entities field
    mesh.entities.push_back(entity);
    return true;
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

    //Can we retrieve info with phys group tag, or is it by chance than physgroupTag = entityTag ?
    // collect the information contained in the gmsh file
    std::vector<std::pair<int, int>> physGroupHandles;
    gmsh::model::getPhysicalGroups(physGroupHandles, 2);

    std::vector<std::pair<int, int>> BCHandles;
    gmsh::model::getPhysicalGroups(BCHandles, 1);

    for(unsigned int i = 0 ; i < BCHandles.size() ; ++i)
    {
        std::string name;
        gmsh::model::getPhysicalName(1, BCHandles[i].second, name);
        std::vector<int> nodesTags;
        std::vector<double> coord;
        gmsh::model::mesh::getNodesForPhysicalGroup(1, BCHandles[i].second,
                                                    nodesTags, coord);
        mesh.nodesTagBoundary[name] = nodesTags;
        //mesh.coordNodesBoundary[name] = coord;
    }

    //Assume 1 2D physical group, might change later
    //New structure ?
    std::vector<int> entitiesTag;
    gmsh::model::getEntitiesForPhysicalGroup(physGroupHandles[0].first,
                                             physGroupHandles[0].second,
                                             entitiesTag);

    //Modify addEntity ?
    std::vector<std::pair<int, int>> entityHandles;
    for(auto entityTag : entitiesTag)
    {
        entityHandles.push_back(std::pair<int, int>(physGroupHandles[0].first, entityTag));
    }

    unsigned int currentOffset = 0;

    for(auto entityHandle : entityHandles)
    {
        if(!addEntity(mesh, entityHandle, currentOffset, intScheme, basisFuncType))
            return false;
    }

    loadNodeData(mesh);

    gmsh::finalize();

    return true;
}
