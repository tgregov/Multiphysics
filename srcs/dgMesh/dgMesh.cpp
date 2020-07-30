#include "dgMesh.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <gmsh.h>

dgMesh::dgMesh(std::string integrationType, std::string basisFuncType):
m_dimension(-1),
m_integrationType(integrationType),
m_basisFuncType(basisFuncType)
{

}

void dgMesh::loadFromFile(std::string fileName)
{
    /**
        NB: In this work, we will only support mesh with a single element type.

        A geometry is represented by a set of geometrical entities. There is normally entities
        for every geometrical object of any dimension in the mesh (points, curves, surfaces and volumes).
        Any entity should belong to **at least** one physical group. Every entity is then discretized into
        a set of elements (which are called elements regardless of the dimension).

        NB: two entities can have the same tag if they are of different dimension, but not physical groups.

        Multiple meshed of different element types can composed an entity, but we will not support them.

        Despite every entity of any dimension being meshed using elements, one can retrieve the "edges"
        and "faces" of an element.
    **/

    //Dummy stream to check if the file exist, because gmsh does not do it.
    std::ifstream dummyStream(fileName);
    if(!dummyStream.is_open())
        throw std::runtime_error("could not open file " + fileName + ".");

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName);

    m_dimension = gmsh::model::getDimension();

    loadPhysicalGroupsAndEntities();

    loadElementsProperty();

    loadElements();

    gmsh::finalize();
}

void dgMesh::displayToConsole() const noexcept
{
    std::cout << m_dimension << "D physical groups: " << std::endl;
    for(PhysicalGroup pg : m_physicalGroupHD)
    {
        std::cout << "\t * " << pg.name << ": \n";
        std::cout << "\t\t - tag: " << pg.tag << "\n";
        std::cout << "\t\t - entities tag: ";
        for(Entity* pEntity : pg.pEntities)
            std::cout << pEntity->mainTag << ", ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << m_dimension - 1 << "D physical groups: " << std::endl;
    for(PhysicalGroup pg : m_physicalGroupLD)
    {
        std::cout << "\t * " << pg.name << ": \n";
        std::cout << "\t\t - tag: " << pg.tag << "\n";
        std::cout << "\t\t - entities tag: ";
        for(Entity* pEntity : pg.pEntities)
            std::cout << pEntity->mainTag << ", ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << m_dimension << "D entities: " << std::endl;
    for(Entity entity : m_entitiesHD)
    {
        std::cout << "\t * " << entity.mainTag << ": \n";
        std::cout << "\t\t - element type: " << entity.pElementProperty->type << "\n";
        std::cout << "\t\t - element name: " << entity.pElementProperty->name << "\n";
    }
    std::cout << std::endl;

    std::cout << m_dimension - 1 << "D entities: " << std::endl;
    for(Entity entity : m_entitiesLD)
    {
        std::cout << "\t * " << entity.mainTag << ": \n";
        std::cout << "\t\t - element type: " << entity.pElementProperty->type << "\n";
        std::cout << "\t\t - element name: " << entity.pElementProperty->name << "\n";
    }
    std::cout << std::endl;
}

void dgMesh::loadPhysicalGroupsAndEntities()
{
    /** We firstly retrieve the physical groups of dimension HD and LD **/

    gmsh::vectorpair physicalGroupsHD;
    gmsh::model::getPhysicalGroups(physicalGroupsHD, m_dimension);
    if(physicalGroupsHD.size() == 0)
    {
        throw std::runtime_error("it seems your mesh does not have any physical groups set for dimension "
                                 + std::to_string(m_dimension) + ".");
    }

    for(std::pair pgHD : physicalGroupsHD)
    {
        std::string pgName;
        gmsh::model::getPhysicalName(pgHD.first, pgHD.second, pgName);

        if(pgName.empty())
        {
            throw std::runtime_error("one of your " + std::to_string(m_dimension) +
                                     "D physical group does not have a name.");
        }

        m_physicalGroupHD.push_back({pgHD.second, pgName, {}});
    }

    gmsh::vectorpair physicalGroupsLD;
    gmsh::model::getPhysicalGroups(physicalGroupsLD, m_dimension - 1);
    if(physicalGroupsLD.size() == 0)
    {
        throw std::runtime_error("it seems your mesh does not have any physical groups set for dimension "
                                 + std::to_string(m_dimension - 1) + ".");
    }

    for(std::pair pgLD : physicalGroupsLD)
    {
        std::string pgName;
        gmsh::model::getPhysicalName(pgLD.first, pgLD.second, pgName);

        if(pgName.empty())
        {
            throw std::runtime_error("one of your " + std::to_string(m_dimension - 1) +
                                     "D physical group does not have a name.");
        }

        m_physicalGroupLD.push_back({pgLD.second, pgName, {}});
    }

    /**
        We then retrieves the entities of dimension HD and LD. Their names are not
        retrieve since we will not use them
    **/

    gmsh::vectorpair entitiesHD;
    gmsh::model::getEntities(entitiesHD, m_dimension);
    if(entitiesHD.size() == 0)
    {
        throw std::runtime_error("it seems your mesh does not have any geometry set for dimension "
                                 + std::to_string(m_dimension) + ".");
    }

    for(std::pair entityHD : entitiesHD)
    {
        m_entitiesHD.push_back({
            entityHD.second,
            -1,
            {},
            nullptr,
            {},
        });
    }

    gmsh::vectorpair entitiesLD;
    gmsh::model::getEntities(entitiesLD, m_dimension -1);
    if(entitiesLD.size() == 0)
    {
        throw std::runtime_error("it seems your mesh does not have any geometry set for dimension "
                                 + std::to_string(m_dimension - 1) + ".");
    }

    for(std::pair entityLD : entitiesLD)
    {
        m_entitiesLD.push_back({
            entityLD.second,
            -1,
            {},
            nullptr,
            {},
        });
    }

    /** We now create the vector of pointers to entities in each physical group **/

    for(PhysicalGroup& pgHD : m_physicalGroupHD)
    {
        std::vector<int> entitiesTagInPG;
        gmsh::model::getEntitiesForPhysicalGroup(m_dimension, pgHD.tag, entitiesTagInPG);

        for(int entityTagInPG : entitiesTagInPG)
        {
            auto it = std::find_if(m_entitiesHD.begin(), m_entitiesHD.end(),
            [entityTagInPG](const Entity& knownEntity)
            {
                return (entityTagInPG == knownEntity.mainTag);
            });

            assert(it != m_entitiesHD.cend());

            pgHD.pEntities.push_back(&*it);
        }
    }

    for(PhysicalGroup& pgLD : m_physicalGroupLD)
    {
        std::vector<int> entitiesTagInPG;
        gmsh::model::getEntitiesForPhysicalGroup(m_dimension - 1, pgLD.tag, entitiesTagInPG);

        for(int entityTagInPG : entitiesTagInPG)
        {
            auto it = std::find_if(m_entitiesLD.begin(), m_entitiesLD.end(),
            [entityTagInPG](const Entity& knownEntity)
            {
                return (entityTagInPG == knownEntity.mainTag);
            });

            assert(it != m_entitiesLD.cend());

            pgLD.pEntities.push_back(&*it);
        }
    }

    /** We now create the vector of pointers to physical group in each entity **/

    for(Entity& entityHD : m_entitiesHD)
    {
        std::vector<int> physicalGroupsTagInEntity;
        gmsh::model::getPhysicalGroupsForEntity(m_dimension, entityHD.mainTag, physicalGroupsTagInEntity);

        for(int pgTagInEntity : physicalGroupsTagInEntity)
        {
            auto it = std::find_if(m_physicalGroupHD.begin(), m_physicalGroupHD.end(),
            [pgTagInEntity](const PhysicalGroup& knownPG)
            {
                return (pgTagInEntity == knownPG.tag);
            });

            assert(it != m_physicalGroupHD.cend());

            entityHD.pPhysicalGroups.push_back(&*it);
        }
    }

    for(Entity& entityLD : m_entitiesLD)
    {
        std::vector<int> physicalGroupsTagInEntity;
        gmsh::model::getPhysicalGroupsForEntity(m_dimension - 1, entityLD.mainTag, physicalGroupsTagInEntity);

        for(int pgTagInEntity : physicalGroupsTagInEntity)
        {
            auto it = std::find_if(m_physicalGroupLD.begin(), m_physicalGroupLD.end(),
            [pgTagInEntity](const PhysicalGroup& knownPG)
            {
                return (pgTagInEntity == knownPG.tag);
            });

            assert(it != m_physicalGroupLD.cend());

            entityLD.pPhysicalGroups.push_back(&*it);
        }
    }
}

void dgMesh::loadElementsProperty()
{
    /** We retrieve all the element property inside the whole mesh **/

    std::vector<int> elementsType;
    gmsh::model::mesh::getElementTypes(elementsType);

    for(int elementType : elementsType)
    {
        ElementProperty elmProperty;
        elmProperty.type = elementType;
        int dimension;
        gmsh::model::mesh::getElementProperties(elementType,
                                                elmProperty.name, dimension,
                                                elmProperty.order, elmProperty.numNodes,
                                                elmProperty.localNodeCoord, elmProperty.numPrimaryNodes);

        gmsh::model::mesh::getIntegrationPoints(elementType, m_integrationType,
                                                elmProperty.intPointsCoord,
                                                elmProperty.intPointsWeigth);

        int dummyNumOrientation;
        gmsh::model::mesh::getBasisFunctions(elementType, elmProperty.intPointsCoord,
                                             m_basisFuncType, elmProperty.basisFuncNumComp,
                                             elmProperty.basisFunctions, dummyNumOrientation);

        assert(dummyNumOrientation == 1);

        gmsh::model::mesh::getBasisFunctions(elementType, elmProperty.intPointsCoord,
                                             "Grad" + m_basisFuncType, elmProperty.basisFuncNumCompGrad,
                                             elmProperty.basisFunctionsGrad, dummyNumOrientation);

        assert(dummyNumOrientation == 1);

        if(dimension == m_dimension)
            m_elementsPropertyHD.push_back(std::move(elmProperty));
        else if(dimension == m_dimension - 1)
            m_elementsPropertyLD.push_back(std::move(elmProperty));
    }

    /** We now determine which element type is in which entity (no hybrid meshes) **/

    for(Entity& entity : m_entitiesHD)
    {
        std::vector<int> elementType;
        gmsh::model::mesh::getElementTypes(elementType, m_dimension, entity.mainTag);

        if(elementType.size() != 1)
            throw std::runtime_error("Hybrid Meshes are currently not handled");

        auto it = std::find_if(m_elementsPropertyHD.begin(), m_elementsPropertyHD.end(),
        [elementType](const ElementProperty& knownElmProp)
        {
            return (elementType[0] == knownElmProp.type);
        });

        assert(it != m_elementsPropertyHD.cend());

        entity.pElementProperty = &*it;
    }

    for(Entity& entity : m_entitiesLD)
    {
        std::vector<int> elementType;
        gmsh::model::mesh::getElementTypes(elementType, m_dimension - 1, entity.mainTag);

        if(elementType.size() != 1)
            throw std::runtime_error("Hybrid Meshes are currently not handled");

        auto it = std::find_if(m_elementsPropertyLD.begin(), m_elementsPropertyLD.end(),
        [elementType](const ElementProperty& knownElmProp)
        {
            return (elementType[0] == knownElmProp.type);
        });

        assert(it != m_elementsPropertyLD.cend());

        entity.pElementProperty = &*it;
    }
}

void dgMesh::loadElements()
{
    /** Element in HD entities should appear inly once **/

    for(Entity& entity : m_entitiesHD)
    {
        entity.subTag = gmsh::model::addDiscreteEntity(m_dimension - 1);
        std::vector<std::size_t> nodesTagPerEdge;
        gmsh::model::mesh::getElementEdgeNodes(entity.pElementProperty->type, nodesTagPerEdge, entity.mainTag);
        int eleTypeSubEntity;
        switch(m_dimension)
        {
            case 1:
                eleTypeSubEntity = gmsh::model::mesh::getElementType("point", entity.pElementProperty->order);
                break;

            case 2:
                eleTypeSubEntity = gmsh::model::mesh::getElementType("line", entity.pElementProperty->order);
                break;

            case 3:
                eleTypeSubEntity = gmsh::model::mesh::getElementType("line", entity.pElementProperty->order);
                break;
        }

        for(auto& elmProp : m_elementsPropertyLD)
        {
            if(elmProp.type == eleTypeSubEntity)
            {
                entity.pFaceEdgePropety = &elmProp;
                break;
            }
        }

        assert(entity.pFaceEdgePropety != nullptr);

        gmsh::model::mesh::addElementsByType(entity.subTag, eleTypeSubEntity, {}, nodesTagPerEdge);

        std::vector<std::size_t> elementTags;
        std::vector<std::size_t> nodeTags;
        gmsh::model::mesh::getElementsByType(entity.pElementProperty->type,
                                             elementTags, nodeTags, entity.mainTag);

        std::vector<double> jacobiansHD;
        std::vector<double> determinantsHD;
        std::vector<double> dummyCoordHD;
        gmsh::model::mesh::getJacobians(entity.pElementProperty->type, entity.pElementProperty->intPointsCoord,
                                        jacobiansHD, determinantsHD, dummyCoordHD, entity.mainTag);

        std::vector<double> jacobiansLD;
        std::vector<double> determinantsLD;
        std::vector<double> dummyCoordLD;
        gmsh::model::mesh::getJacobians(entity.pFaceEdgePropety->type, entity.pFaceEdgePropety->intPointsCoord,
                                        jacobiansLD, determinantsLD, dummyCoordLD, entity.subTag);

        std::vector<double> baryCenters;
        gmsh::model::mesh::getBarycenters(entity.pElementProperty->type, entity.mainTag, false, true, baryCenters);

        unsigned int nNodesElm = entity.pElementProperty->numNodes;
        unsigned int nElements = elementTags.size();
        unsigned int nNodesFaceEdge = entity.pFaceEdgePropety->numNodes;
        unsigned int nFaceEdges = nodesTagPerEdge.size()/nNodesFaceEdge;
        unsigned int nFaceEdgesPerElm = nNodesElm;
        unsigned int nGPHD = entity.pElementProperty->intPointsWeigth.size();
        unsigned int nGPLD = entity.pFaceEdgePropety->intPointsWeigth.size();

        for(std::size_t e = 0 ; e < nElements ; ++e)
        {
            Element elm = {};
            elm.tag         = elementTags[e];
            elm.nodesTag    = std::vector<std::size_t>(nodeTags.begin() + nNodesElm*e,
                                                       nodeTags.begin() + nNodesElm*(e + 1));
            elm.determinant = std::vector<double>(determinantsHD.begin() + nGPHD*e,
                                                  determinantsHD.begin() + nGPHD*(e + 1));
            elm.jacobian    = std::vector<double>(jacobiansHD.begin() + 9*nGPHD*e,
                                                  jacobiansHD.begin() + 9*nGPHD*(e + 1));
            elm.pEntity     = &entity;

            for(std::size_t nodeTag : elm.nodesTag)
            {
                std::vector<double> coord, dummyParametricCoord;
                gmsh::model::mesh::getNode(nodeTag, coord, dummyParametricCoord);
                elm.nodesCoord.push_back(coord);
            }

            std::vector<double> elementBarycenters(baryCenters.begin() + 3*e,
                                                   baryCenters.begin() + 3*(e + 1));

            for(unsigned int i = 0 ; i < nFaceEdgesPerElm ; ++i)
            {
                FaceEdge faceEdge = {};

                faceEdge.nodesTag = std::vector<std::size_t>(nodesTagPerEdge.begin() + nNodesFaceEdge*nFaceEdgesPerElm*e + nNodesFaceEdge*i,
                                                             nodesTagPerEdge.begin() + nNodesFaceEdge*nFaceEdgesPerElm*e + nNodesFaceEdge*(i + 1));

                for(std::size_t nodeTag : faceEdge.nodesTag)
                {
                    std::vector<double> coord, dummyParametricCoord;
                    gmsh::model::mesh::getNode(nodeTag, coord, dummyParametricCoord);
                    faceEdge.nodesCoord.push_back(coord);
                }

                faceEdge.determinant = std::vector<double>(determinantsLD.begin() + nFaceEdgesPerElm*nGPLD*e + nGPLD*i,
                                                           determinantsLD.begin() + nFaceEdgesPerElm*nGPLD*e + nGPLD*(i + 1));

                faceEdge.jacobian = std::vector<double>(jacobiansLD.begin() + nFaceEdgesPerElm*9*nGPLD*e + 9*nGPLD*i,
                                                        jacobiansLD.begin() + nFaceEdgesPerElm*9*nGPLD*e + 9*nGPLD*(i + 1));

                computeFaceEdgeNormal(faceEdge, elementBarycenters);

                m_faceEdges.push_back(std::move(faceEdge));
                elm.pFaceEdges.push_back(&m_faceEdges.back());
            }

            m_elementsHD.push_back(std::move(elm));
            entity.pElements.push_back(&m_elementsHD.back());

            Element& finalElm = m_elementsHD.back();
            for(FaceEdge* pFaceEdge : finalElm.pFaceEdges)
            {
                pFaceEdge->parentElementHD = &finalElm;
            }
        }
    }
}

void dgMesh::computeFaceEdgeNormal(FaceEdge& faceEdge, const std::vector<double>& elementBarycenter)
{
    std::vector<double> normal;

    //TO DO: 3D version and 1D version
//    double edgeLength = std::sqrt((faceEdge.nodesCoord[0][0] - faceEdge.nodesCoord[1][0])
//                                  *(faceEdge.nodesCoord[0][0] - faceEdge.nodesCoord[1][0])
//                                 + (faceEdge.nodesCoord[0][1] - faceEdge.nodesCoord[1][1])
//                                  *(faceEdge.nodesCoord[0][1] - faceEdge.nodesCoord[1][1])
//                                 + (faceEdge.nodesCoord[0][2] - faceEdge.nodesCoord[1][2])
//                                  *(faceEdge.nodesCoord[0][2] - faceEdge.nodesCoord[1][2]));

    // compute the normal
    // if A:(x1, y1) and B:(x2, y2), then AB = (x2 - x1, y2 - y1) and a
    // normal is given by n = (y2 - y1, x1 - x2)
    double nx = faceEdge.nodesCoord[1][1] - faceEdge.nodesCoord[0][1];
    double ny = faceEdge.nodesCoord[0][0] - faceEdge.nodesCoord[1][0];
    double norm = sqrt(ny*ny + nx*nx);

    // unfortunately, nodes per edge in nodes vector are not always in
    // the same order (clockwise vs anticlockwise) => we need to check
    // the orientation
    double vx = elementBarycenter[0] - (faceEdge.nodesCoord[1][0] + faceEdge.nodesCoord[0][0])/2;
    double vy = elementBarycenter[1] - (faceEdge.nodesCoord[1][1] + faceEdge.nodesCoord[0][1])/2;

    if(nx*vx + ny*vy > 0)
    {
        nx = -nx;
        ny = -ny;
    }

    // normalize the normal components
    normal.push_back(nx/norm);
    normal.push_back(ny/norm);

    faceEdge.normal = std::move(normal);
}
