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
    gmsh::open(fileName);

    m_dimension = gmsh::model::getDimension();

    loadPhysicalGroupsAndEntities();

    loadElementsProperty();

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
            std::cout << pEntity->tag << ", ";
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
            std::cout << pEntity->tag << ", ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << m_dimension << "D entities: " << std::endl;
    for(Entity entity : m_entitiesHD)
    {
        std::cout << "\t * " << entity.tag << ": \n";
        std::cout << "\t\t - element type: " << entity.pElementProperty->type << "\n";
        std::cout << "\t\t - element name: " << entity.pElementProperty->name << "\n";
    }
    std::cout << std::endl;

    std::cout << m_dimension - 1 << "D entities: " << std::endl;
    for(Entity entity : m_entitiesLD)
    {
        std::cout << "\t * " << entity.tag << ": \n";
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
                return (entityTagInPG == knownEntity.tag);
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
                return (entityTagInPG == knownEntity.tag);
            });

            assert(it != m_entitiesLD.cend());

            pgLD.pEntities.push_back(&*it);
        }
    }

    /** We now create the vector of pointers to physical group in each entity **/

    for(Entity& entityHD : m_entitiesHD)
    {
        std::vector<int> physicalGroupsTagInEntity;
        gmsh::model::getPhysicalGroupsForEntity(m_dimension, entityHD.tag, physicalGroupsTagInEntity);

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
        gmsh::model::getPhysicalGroupsForEntity(m_dimension - 1, entityLD.tag, physicalGroupsTagInEntity);

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
        gmsh::model::mesh::getElementTypes(elementType, m_dimension, entity.tag);

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
        gmsh::model::mesh::getElementTypes(elementType, m_dimension - 1, entity.tag);

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
