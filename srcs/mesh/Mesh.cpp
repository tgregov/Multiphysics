/**
 * \file Mesh.cpp
 * \brief Implementation of the required function to load a Mesh structure from file.
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <gmsh.h>
#include "Mesh.hpp"
#include "../utils/utils.hpp"


static void loadNodeData(Mesh& mesh)
{

    unsigned int numNodes = 0;
    std::vector<std::size_t> elementTags;
    std::vector<unsigned int> elementNumNodes;
    std::vector<std::size_t> nodeTags;
    std::vector<std::vector<double>> coord;

    // loop over the elements
    for(size_t elm = 0 ; elm < mesh.elements.size() ; ++elm)
    {
        elementTags.push_back(mesh.elements[elm].elementTag);
        elementNumNodes.push_back(mesh.elements[elm].nodeTags.size());

        for(size_t n = 0 ; n < mesh.elements[elm].nodeTags.size() ; ++n)
        {
            std::vector<double> temp;
            temp.push_back(mesh.elements[elm].nodesCoord[n][0]);
            temp.push_back(mesh.elements[elm].nodesCoord[n][1]);
            coord.push_back(temp);

            nodeTags.push_back(mesh.elements[elm].nodeTags[n]);
            numNodes++;

        }
    }

    mesh.nodeData.numNodes = numNodes;
    mesh.nodeData.elementTags = elementTags;
    mesh.nodeData.elementNumNodes = elementNumNodes;
    mesh.nodeData.nodeTags = nodeTags;
    mesh.nodeData.coord = coord;
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

            gmsh::model::mesh::getIntegrationPoints(eleTypes[i], intScheme,
                                                    elementProperty.intPoints,
                                                    elementProperty.intWeigths);

            gmsh::model::mesh::getBasisFunctions(eleTypes[i],
                                                 elementProperty.intPoints,
                                                 basisFuncType,
                                                 elementProperty.numComp,
                                                 elementProperty.basisFunc);

            int dummyNumComp;
            gmsh::model::mesh::getBasisFunctions(eleTypes[i],
                                                 elementProperty.intPoints,
                                                 std::string("Grad" + basisFuncType),
                                                 dummyNumComp,
                                                 elementProperty.basisFuncGrad);


            // add the products prodFunc[k][i,j] = w_k*l_i(x_k)*l_j(x_k)
            // add the products pondFunc[k][i] = w_k*l_i(x_k)
            elementProperty.nGP = elementProperty.intPoints.size()/3;
            elementProperty.nSF = elementProperty.basisFunc.size()
                                    /elementProperty.nGP;
            for(unsigned int k = 0 ; k < elementProperty.nGP ; ++k)
            {
                std::vector<double> wll;
                std::vector<double> wl;

                for(unsigned int i = 0 ; i < elementProperty.nSF ; ++i)
                {

                    wl.push_back(elementProperty.intWeigths[k]
                            *elementProperty.basisFunc[elementProperty.nSF*k + i]);

                    for(unsigned int j = i ; j < elementProperty.nSF ; ++j)
                    {
                        if(k == 0)
                        {
                            elementProperty.IJ.push_back(
                                std::pair<unsigned int, unsigned int>(i, j));
                        }

                        wll.push_back(elementProperty.intWeigths[k]
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

            for(unsigned int l = 0 ;
                    l < elementProperty.nSF*(elementProperty.nSF+1)/2 ; ++l)
            {
                double  sum = 0.0;
                for (unsigned int k = 0 ; k < elementProperty.nGP ; ++k)
                {
                    sum += elementProperty.prodFunc[k][l];
                }

                lalb[elementProperty.IJ[l].first][elementProperty.IJ[l].second] = sum;
                if(elementProperty.IJ[l].first != elementProperty.IJ[l].second)
                {
                    lalb[elementProperty.IJ[l].second][elementProperty.IJ[l].first]
                        = sum;
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
 * \param nodesTagsEdge Nodes tag per edge of the element.
 * \param determinantLD Determinant associated with the edge of the element.
 * \param nodesTagsElement Nodes tag of the element.
 */
static void addEdge(Element& element, std::vector<int> nodesTagsEdge,
                    std::vector<double> determinantLD,
                    const std::vector<int>& nodesTagsElement)
{

    Edge edge;
    // should be before movement !
    for(unsigned int i = 0 ; i < nodesTagsEdge.size() ; ++i)
    {
        for(unsigned int j = 0 ; j < nodesTagsElement.size() ; ++j)
        {
            if(nodesTagsEdge[i] == nodesTagsElement[j]) // we should always find it!
                edge.offsetInElm.push_back(j);

        }
    }
    edge.nodeTags = std::move(nodesTagsEdge);
    edge.determinantLD = std::move(determinantLD);

    element.edges.push_back(edge);
    // we still need to add its tag
}


/**
 * \brief Compute the outward normal of an edge.
 * \param edge The edge of which the normal is computed.
 * \param meshDim Dimension of the mesh (1 or 2)
 * \param baryCenter  Barycenter of the parent element.
 */
static void computeEdgeNormalCoord(Edge& edge, unsigned int meshDim,
                              const std::vector<double>& baryCenter)
{
    std::vector<double> normal;

    std::vector<double> coord1, coord2;
    for(unsigned int i = 0 ; i < edge.nodeTags.size() ; ++i)
    {
            std::vector<double> coord, dummyParametricCoord;
            gmsh::model::mesh::getNode(edge.nodeTags[i], coord,
                                        dummyParametricCoord);

            edge.nodeCoordinate.push_back(coord);

            // normally first and second node tags are vertex
            if(i == 0)
                coord1 = std::move(coord);
            else if(i == 1)
                coord2 = std::move(coord);
    }

    edge.length=sqrt((coord1[0]-coord2[0])*(coord1[0]-coord2[0])
                             +(coord1[1]-coord2[1])*(coord1[1]-coord2[1])
                             +(coord1[2]-coord2[2])*(coord1[2]-coord2[2]));

    switch(meshDim)
    {
        case 1:
            // normals are either 1 or -1
            if(coord1[0]-baryCenter[0] < 0)
                normal.push_back(-1);
            else
                normal.push_back(1);
                edge.normal = std::move(normal);
                break;

        case 2:
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
            normal.push_back(nx/norm);
            normal.push_back(ny/norm);

            edge.normal = std::move(normal);
            break;
    }
}

/**
 * \brief Find the second position of the edge in the entity.
 * \param entity The parent entity.
 * \param currentEdge The edge of which you want to find the neighbor.
 * \param edgePos Index of the edge in its element.
 */
static void findInFrontEdge(Mesh& mesh, Edge& currentEdge, unsigned int edgePos)
{
    // try to find the edge in another element. If an edge is not a boundary, then
    // it has necessarily an "edge in front", which will be found
    // (the check is done an all previously computed elements).
    unsigned int elVecSize = mesh.elements.size();
    unsigned int nEdgePerEl = mesh.elements[0].edges.size();

    #pragma omp parallel default(none) shared(mesh, elVecSize, nEdgePerEl, currentEdge, edgePos)
    {
        #pragma omp for
        for(unsigned int elm = 0 ; elm < elVecSize ; ++elm)
        {
            for(unsigned int k = 0 ; k < nEdgePerEl ; ++k)
            {
                std::vector<int> nodesTagsInfront
                    = mesh.elements[elm].edges[k].nodeTags;
                std::vector<unsigned int> permutation1, permutation2;

                // if the nodes tags of an edge is the permutation of the nodes tags
                //  of another, we have found the "edge in front"
                if(isPermutation(currentEdge.nodeTags, nodesTagsInfront,
                                    permutation1, permutation2))
                {
                    #pragma omp critical
                    {
                        currentEdge.edgeInFront
                            =  std::pair<unsigned int, unsigned int>(elm, k);
                        mesh.elements[elm].edges[k].edgeInFront
                            = std::pair<unsigned int, unsigned int>
                            (elVecSize, edgePos);
                        currentEdge.nodeIndexEdgeInFront = std::move(permutation1);
                        mesh.elements[elm].edges[k].nodeIndexEdgeInFront
                            = std::move(permutation2);
                    }
                    #pragma omp cancel for
                }
            }
            #pragma omp cancellation point for
        }
    }
}

/**
 * \brief Check if an edge lies on a boundary (and which one)
 * \param nodesTagBoundaries Map which stores the tags of the nodes belonging to a
 *  certain boundary.
 * \param edge The edge which we check if it is a boundary.
 */
static bool IsBoundary(const std::map<std::string,
                            std::vector<std::size_t>>& nodesTagBoundaries, Edge& edge)
{
    // check if an edge is on the boundary of the domain.
    for(std::pair<std::string, std::vector<std::size_t>> nodeTagBoundary : nodesTagBoundaries)
    {
        for(unsigned int j = 0 ; j < edge.nodeTags.size() ; ++j)
        {
            // if one of the nodes tags of the edge is not found inside the boundary
            // then the edge does not belong to that boundary
            if(!std::count(nodeTagBoundary.second.begin(),
                            nodeTagBoundary.second.end(), edge.nodeTags[j]))
                break;

            // ok the edge belongs to that boundary.
            if(j == edge.nodeTags.size() - 1)
            {
                edge.bcName = nodeTagBoundary.first;
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
 * \param eleTypeHD Type of the element.
 * \param eleTypeLD Type of the element's edges.
 * \param jacobiansHD Jacobian matrix of the variable change of the element
 * evaluated at each Gauss points.
 * \param determinantsHD Determinant of the variable change of the element
 * evaluated at each Gauss points.
 * \param determinantsLD Determinant of the variable change of the element's edges
 * evaluated at each Gauss points.
 * \param nGPLD Number of Gauss point for integration.
 * \param offsetInU Offset of the element in the u unknown vector.
 * \param nodesTagsPerEdge Node tags of the element, per edge.
 * \param nodesTags Node tags of the element.
 * \param elementBarycenter Element barycenter coordinate.
 */
static void addElement(Mesh& mesh, int elementTag, int eleTypeHD,
                        int eleTypeLD, std::vector<double> jacobiansHD,
                        std::vector<double> determinantsHD,
                        std::vector<double> determinantsLD,
                        unsigned int nGPLD, unsigned int offsetInU,
                        std::vector<int> nodesTagsPerEdge,
                        std::vector<int> nodesTags,
                        const std::vector<double>& elementBarycenter)
{
    // fill an element structure
    Element element;
    element.elementTag = elementTag;
    element.elementTypeHD = eleTypeHD;
    element.elementTypeLD = eleTypeLD;

    element.offsetInU = offsetInU;

    element.determinantHD = std::move(determinantsHD);
    element.jacobianHD = std::move(jacobiansHD);
    element.nodeTags = std::move(nodesTags);

    for(unsigned int i = 0 ; i < element.nodeTags.size() ; ++i)
    {
        std::vector<double> coord, dummyParametricCoord;
        gmsh::model::mesh::getNode(element.nodeTags[i], coord, dummyParametricCoord);

        element.nodesCoord.push_back(coord);
    }

    // compute the number of edge and of nodes per edge of that element
    unsigned int nEdge = determinantsLD.size()/nGPLD;
    unsigned int nNodePerEdge = mesh.elementProperties.at(eleTypeLD).numNodes;

    // for each edge, we add an edge to the element.edges field
    for(unsigned int i = 0 ; i < nEdge ; ++i)
    {
        // get the nodes tags associated with that particular edge
        std::vector<int> nodesTagsEdge(nodesTagsPerEdge.begin() + nNodePerEdge*i,
                                        nodesTagsPerEdge.begin()
                                        + nNodePerEdge*(i + 1));

        // get the determinants associated with that particular edge
        std::vector<double> determinantsEdgeLD(determinantsLD.begin()
                                                + nGPLD*i, determinantsLD.begin()
                                                + nGPLD*(i + 1));

        // add the edge to the element.edges field
        addEdge(element, std::move(nodesTagsEdge), std::move(determinantsEdgeLD),
                    element.nodeTags);
        if(mesh.elements.size() != 0)
        {
            // check if the edge is on the boundary of the domain
            // if not, try to find the edge neighbour in a previously computed element
            if(!IsBoundary(mesh.nodesTagBoundary, element.edges[i]))
                findInFrontEdge(mesh, element.edges[i], i);

        }
        else
        {
            IsBoundary(mesh.nodesTagBoundary, element.edges[i]);
        }

        // compute the normal of the edge and get the nodes coordinates of the edge
        computeEdgeNormalCoord(element.edges[i], mesh.dim, elementBarycenter);

        // compute some essential matrix for the DG-FEM method
        Eigen::SparseMatrix<double> dMs(mesh.elementProperties.at(eleTypeHD).nSF,
                                        mesh.elementProperties.at(eleTypeHD).nSF);
        std::vector<Eigen::Triplet<double>> indices;

        for(size_t nA = 0 ; nA <  element.edges[i].nodeTags.size() ; ++nA)
        {
            for(size_t nB = 0 ; nB <  element.edges[i].nodeTags.size() ; ++nB)
            {
                indices.push_back(Eigen::Triplet<double>
                    (element.edges[i].offsetInElm[nA],
                        element.edges[i].offsetInElm[nB],
                        mesh.elementProperties.at(eleTypeLD).lalb[nA][nB]));
            }
        }
        dMs.setFromTriplets(indices.begin(), indices.end());
        element.dM.push_back(dMs);
    }

    // add the element to the mesh
    mesh.elements.push_back(element);
}


/**
 * \brief Add an entity to a certain 2D mesh (filling the required fields).
 * \param mesh The parent mesh.
 * \param entityTag Entity tag to add.
 * \param currentOffset Offset in u unknown vector.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 * \return true if then entity was added flawlessly, false otherwise.
 */
static bool buildMesh(Mesh& mesh, int entityTag, unsigned int& currentOffset,
                      const std::string& intScheme, const std::string& basisFuncType)
{
    // fill an entity structure
    mesh.entityTagHD = entityTag;

    // add an entity for the mesh.dim-1 D elements created (edges of the elements)
    int c = gmsh::model::addDiscreteEntity(1);
    mesh.entityTagLD = c;

    // get the element types in the current mesh.dim D entity
    std::vector<int> eleTypesHD;
    gmsh::model::mesh::getElementTypes(eleTypesHD, mesh.dim,entityTag);

    loadElementProperties(mesh.elementProperties, eleTypesHD,
                          intScheme, basisFuncType);

    // loop over the element types in the current mesh.dim D entity
    for(auto eleTypeHD : eleTypesHD)
    {
        // get the elements of type eleTypeHD and their nodes tags.
        std::vector<std::size_t> nodeTags, elementTags;
        gmsh::model::mesh::getElementsByType(eleTypeHD, elementTags, nodeTags,
                                             entityTag);

        // get the elements nodes tags "per edge".
        std::vector<std::size_t> nodesTagPerEdge;
        gmsh::model::mesh::getElementEdgeNodes(eleTypeHD, nodesTagPerEdge,
                                                entityTag);

        unsigned int numNodes = mesh.elementProperties[eleTypeHD].numNodes;
        unsigned int order = mesh.elementProperties[eleTypeHD].order;

        if(intScheme == "Lagrange" && order > 7)
        {
            std::cerr << "Lagrange polynomials are not stable with an order"
                      << " superior to 7"
                      << std::endl;

            return false;
        }

        // get the elements barycenters (for normal computation)
        std::vector<double> baryCenters;
        gmsh::model::mesh::getBarycenters(eleTypeHD, entityTag, false, true,
                                            baryCenters);

        // add mesh.dim-1 D entity to store all the lines associated to elements
        // of the same order
        int eleTypeLD;
        switch(mesh.dim)
        {
            case 1:
                eleTypeLD = gmsh::model::mesh::getElementType("point", order);
                break;

            case 2:
                eleTypeLD = gmsh::model::mesh::getElementType("line", order);
                break;
        }

        // creation of the mesh.dim-1 D elements
        gmsh::model::mesh::addElementsByType(c, eleTypeLD, {}, nodesTagPerEdge);

        loadElementProperties(mesh.elementProperties, std::vector<int>(1, eleTypeLD),
                              intScheme, basisFuncType);

        // we then load jacobian matrices and their determinants of the mesh.dim D
        // and mesh.dim -1 D elements (useful for M, Sx, Sy, Sz)
        std::vector<double> jacobiansHD, determinantsHD, dummyPointsHD;
        gmsh::model::mesh::getJacobians(eleTypeHD,
                                        mesh.elementProperties[eleTypeHD].intPoints,
                                        jacobiansHD, determinantsHD, dummyPointsHD,
                                        entityTag);

        std::vector<double> dummyJacobiansLD, determinantsLD, dummyPointsLD;
        gmsh::model::mesh::getJacobians(eleTypeLD,
                                        mesh.elementProperties[eleTypeLD].intPoints,
                                        dummyJacobiansLD, determinantsLD,
                                        dummyPointsLD, c);

        // computation of the number of mesh.dim D elements, number of gauss points
        // for mesh.dim and mesh.dim-1 D elements, the number of edges per mesh.dim D
        // elements and the number of nodes per edge.
        unsigned int nElements = nodeTags.size()
                                /mesh.elementProperties[eleTypeHD].numNodes;
        unsigned int nGPHD = mesh.elementProperties[eleTypeHD].nGP;
        unsigned int nGPLD = mesh.elementProperties[eleTypeLD].nGP;
        unsigned int nEdgePerElement = determinantsLD.size()/(nGPLD*nElements);
        unsigned int nNodesPerEdge = mesh.elementProperties[eleTypeLD].numNodes;

        // loop over each mesh.dim D elements
        unsigned ratio, currentDecade = 0;
        for(unsigned int i = 0 ; i < elementTags.size() ; ++i)
        {

            // display progress
            ratio = int(100*double(i)/double(elementTags.size()));
            if(ratio >= currentDecade)
            {
                std::cout   << "\r" << "Entity [" << mesh.entityTagHD << "]: "
                            << ratio << "% of the elements computed"
                            << std::flush;
                currentDecade = ratio + 1;
            }

            // get jacobians and determinant associated with that particular element
            std::vector<double> jacobiansElementHD(
                                            jacobiansHD.begin() + 9*nGPHD*i,
                                            jacobiansHD.begin() + 9*nGPHD*(1 + i));
            std::vector<double> determinantsElementHD(
                                            determinantsHD.begin() + nGPHD*i,
                                            determinantsHD.begin() + nGPHD*(1 + i));
            std::vector<double> determinantElementLD(
                                            determinantsLD.begin() + nEdgePerElement
                                                                    *nGPLD*i,
                                            determinantsLD.begin() + nEdgePerElement
                                                                    *nGPLD*(1 + i));

            // get nodes tags and nodes tags "per edge" associated with that
            // particular element
            std::vector<int> nodeTagsElement(nodeTags.begin() + numNodes*i,
                                             nodeTags.begin() + numNodes*(1 + i));

            std::vector<int> nodesTagPerEdgeElement(
                                nodesTagPerEdge.begin()
                                + nNodesPerEdge*nEdgePerElement*i,
                                nodesTagPerEdge.begin()
                                + nNodesPerEdge*nEdgePerElement*(i + 1));

            // get the barycenter of that particular element
            std::vector<double> elementBarycenter(baryCenters.begin() + 3*i,
                                                    baryCenters.begin() + 3*(i + 1));

            // offset of the element in the unknown vector
            unsigned int elementOffset = numNodes;

            double dx;

            // add the element to the entity
            addElement(mesh, elementTags[i], eleTypeHD, eleTypeLD,
                        std::move(jacobiansElementHD),
                        std::move(determinantsElementHD),
                        std::move(determinantElementLD),
                        nGPLD, currentOffset,
                        std::move(nodesTagPerEdgeElement),
                        std::move(nodeTagsElement),
                        elementBarycenter);

            currentOffset += elementOffset;
        }

         std::cout  << "\r" << "Entity [" << mesh.entityTagHD << "]: "
                    << "100% of the elements computed" << std::flush << std::endl;
    }

    return true;
}


/**
 * \brief Compute the mesh dimension.
 * \return The mesh dimension (1, 2 or 3).
 */
static unsigned short getMeshDim()
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
        }
    }

    return elementDim;
}

// documentation in .hpp file
std::vector<int> getTags(const Mesh& mesh)
{
    std::vector<int> listTags;

    // loop over the nodes
    for(unsigned int elm = 0 ; elm < mesh.elements.size() ; ++elm)
    {
        for(unsigned int n = 0 ; n < mesh.elements[elm].nodeTags.size() ; ++n)
        {
            listTags.push_back(mesh.elements[elm].nodeTags[n]);
        }
    }

    return listTags;
}

// documentation in .hpp file
bool readMesh(Mesh& mesh, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType)
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    std::ifstream file(fileName);
    if(file.is_open())
        file.close();
    else
    {
        std::cerr << "File: " << fileName << " does not exist!" << std::endl;
        return false;
    }
    gmsh::open(fileName);

    // check that the mesh is not 3D
    mesh.dim = getMeshDim();
    if(mesh.dim == 3)
    {
        std::cerr << "3D meshes unsupported!" << std::endl;
        return false;
    }

    // we retrieve the tags of the physical groups of dimension mesh.dim and
    // mesh.dim-1
    std::vector<std::pair<int, int>> physGroupHandles;
    gmsh::model::getPhysicalGroups(physGroupHandles, mesh.dim);

    std::vector<std::pair<int, int>> BCHandles;
    gmsh::model::getPhysicalGroups(BCHandles, mesh.dim - 1);

    // physical groups of dimension mesh.dim-1 are boundary conditions.
    // we then retrieve and store the tag of each node in those boundaries.
    for(auto BCHandle : BCHandles)
    {
        std::string name;
        gmsh::model::getPhysicalName(mesh.dim - 1, BCHandle.second, name);
        std::vector<std::size_t> nodesTags;
        std::vector<double> dummyCoord;
        gmsh::model::mesh::getNodesForPhysicalGroup(mesh.dim - 1, BCHandle.second,
                                                    nodesTags, dummyCoord);
        mesh.nodesTagBoundary[name] = nodesTags;
    }

    // we assume that a physical group contains only one entity and we retrieve them.
    std::vector<int> entitiesTag;
    for(auto physGroupHandle : physGroupHandles)
    {
        std::vector<int> entityTag;
        gmsh::model::getEntitiesForPhysicalGroup(physGroupHandle.first,
                                             physGroupHandle.second,
                                             entityTag);

        entitiesTag.push_back(entityTag[0]);
    }


    unsigned int currentOffset = 0;

    // we add each identified entity to the mesh.
    if(entitiesTag.size() > 1)
    {
        std::cerr << "Multiple " << mesh.dim << "D entities in " << mesh.dim
                  << "D meshes " << "currently not supported" <<std::endl;

        return false;
    }

    for(auto entityTag : entitiesTag)
    {
        if(!buildMesh(mesh, entityTag, currentOffset, intScheme, basisFuncType))
            return false;
    }

    loadNodeData(mesh);

    gmsh::finalize();

    return true;
}
