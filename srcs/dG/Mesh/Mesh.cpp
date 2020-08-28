#include "Mesh.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <gmsh.h>

namespace dG
{
    Mesh::Mesh(std::string integrationType, std::string basisFuncType):
    m_dimension(-1),
    m_integrationType(integrationType),
    m_basisFuncType(basisFuncType)
    {

    }

    void Mesh::loadFromFile(std::string fileName)
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

        checkIfFaceIsBoundary();

        associateFaces();

        gmsh::finalize();
    }

    void Mesh::displayToConsole() const noexcept
    {
        std::cout << m_dimension << "D physical groups: " << std::endl;
        for(const PhysicalGroup& pg : m_physicalGroupHD)
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
        for(const PhysicalGroup& pg : m_physicalGroupLD)
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
        for(const Entity& entity : m_entitiesHD)
        {
            std::cout << "\t * " << entity.mainTag << ": \n";
            std::cout << "\t\t - element type: " << entity.pElementProperty->type << "\n";
            std::cout << "\t\t - element name: " << entity.pElementProperty->name << "\n";
        }
        std::cout << std::endl;

        std::cout << m_dimension - 1 << "D entities: " << std::endl;
        for(const Entity& entity : m_entitiesLD)
        {
            std::cout << "\t * " << entity.mainTag << ": \n";
            std::cout << "\t\t - element type: " << entity.pElementProperty->type << "\n";
            std::cout << "\t\t - element name: " << entity.pElementProperty->name << "\n";
        }
        std::cout << std::endl;
    }

    void Mesh::loadPhysicalGroupsAndEntities()
    {
        /** We firstly retrieve the physical groups of dimension HD and LD **/

        gmsh::vectorpair physicalGroupsHD;
        gmsh::model::getPhysicalGroups(physicalGroupsHD, m_dimension);
        if(physicalGroupsHD.size() == 0)
        {
            throw std::runtime_error("it seems your mesh does not have any physical groups set for dimension "
                                     + std::to_string(m_dimension) + ".");
        }

        m_physicalGroupHD.resize(physicalGroupsHD.size());
        for(std::size_t i = 0 ; i < m_physicalGroupHD.size() ; ++i)
        {
            std::pair pgHD = physicalGroupsHD[i];

            std::string pgName;
            gmsh::model::getPhysicalName(pgHD.first, pgHD.second, pgName);

            if(pgName.empty())
            {
                throw std::runtime_error("one of your " + std::to_string(m_dimension) +
                                         "D physical group does not have a name.");
            }

            PhysicalGroup physGroup;
            physGroup.tag = pgHD.second;
            physGroup.name = pgName;

            m_physicalGroupHD[i] = std::move(physGroup);
        }

        gmsh::vectorpair physicalGroupsLD;
        gmsh::model::getPhysicalGroups(physicalGroupsLD, m_dimension - 1);
        if(physicalGroupsLD.size() == 0)
        {
            throw std::runtime_error("it seems your mesh does not have any physical groups set for dimension "
                                     + std::to_string(m_dimension - 1) + ".");
        }

        m_physicalGroupLD.resize(physicalGroupsLD.size());
        for(std::size_t i = 0 ; i < m_physicalGroupLD.size() ; ++i)
        {
            std::pair pgLD = physicalGroupsLD[i];

            std::string pgName;
            gmsh::model::getPhysicalName(pgLD.first, pgLD.second, pgName);

            if(pgName.empty())
            {
                throw std::runtime_error("one of your " + std::to_string(m_dimension - 1) +
                                         "D physical group does not have a name.");
            }

            PhysicalGroup physGroup;
            physGroup.tag = pgLD.second;
            physGroup.name = pgName;

            m_physicalGroupLD[i] = std::move(physGroup);
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

        m_entitiesHD.resize(entitiesHD.size());
        for(std::size_t i = 0 ; i < m_entitiesHD.size() ; ++i)
        {
            std::pair entityHD = entitiesHD[i];

            Entity entity;
            entity.mainTag = entityHD.second;

            m_entitiesHD[i] = std::move(entity);
        }

        gmsh::vectorpair entitiesLD;
        gmsh::model::getEntities(entitiesLD, m_dimension -1);
        if(entitiesLD.size() == 0)
        {
            throw std::runtime_error("it seems your mesh does not have any geometry set for dimension "
                                     + std::to_string(m_dimension - 1) + ".");
        }

        m_entitiesLD.resize(entitiesLD.size());
        for(std::size_t i = 0 ; i < m_entitiesLD.size() ; ++i)
        {
            std::pair entityLD = entitiesLD[i];

            Entity entity;
            entity.mainTag = entityLD.second;

            m_entitiesLD[i] = std::move(entity);
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

    void Mesh::loadElementsProperty()
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

    void Mesh::loadElements()
    {
        /**
            Get all the elementTags across all entities (assume one elementType per entities!)
            and resize the vectors of elements LD.
        **/

        std::vector<int> dummyElementTypeLD;
        std::vector<std::vector<std::size_t>> dummyElementTagsLD;
        std::vector<std::vector<std::size_t>> dummyNodesTagsLD;
        gmsh::model::mesh::getElements(dummyElementTypeLD, dummyElementTagsLD, dummyNodesTagsLD, m_dimension - 1);

        std::size_t totalElementLDNumber = 0;
        for(auto elmTags : dummyElementTagsLD)
            totalElementLDNumber += elmTags.size();

        m_elementsLD.resize(totalElementLDNumber);

        /**
            Fill the elements LD vector with th required infos (we do not load determinant
            and jacobian as it will be loaded in the required faceEdge.
        **/

        std::size_t elmLDCounter = 0;

        for(Entity& entity : m_entitiesLD)
        {
            std::vector<std::size_t> elementTags;
            std::vector<std::size_t> nodeTags;
            gmsh::model::mesh::getElementsByType(entity.pElementProperty->type,
                                                 elementTags, nodeTags, entity.mainTag);


            unsigned int nNodesElm = entity.pElementProperty->numNodes;
            unsigned int nElements = elementTags.size();


            for(std::size_t e = 0 ; e < nElements ; ++e)
            {
                Element elm = {};
                elm.tag         = elementTags[e];
                elm.nodesTag    = std::vector<std::size_t>(nodeTags.begin() + nNodesElm*e,
                                                           nodeTags.begin() + nNodesElm*(e + 1));
                elm.pEntity     = &entity;

                for(std::size_t nodeTag : elm.nodesTag)
                {
                    std::vector<double> coord, dummyParametricCoord;
                    gmsh::model::mesh::getNode(nodeTag, coord, dummyParametricCoord);
                    elm.nodesCoord.push_back(coord);
                }

                m_elementsLD[elmLDCounter] = std::move(elm);
                entity.pElements.push_back(&m_elementsLD.back());

                elmLDCounter++;
            }
        }

        /**
            Get all the elementTags across all entities (assume one elementType per entities!)
            and resize the vectors of elements HD.
        **/

        std::vector<int> dummyElementType;
        std::vector<std::vector<std::size_t>> dummyElementTags;
        std::vector<std::vector<std::size_t>> dummyNodesTags;
        gmsh::model::mesh::getElements(dummyElementType, dummyElementTags, dummyNodesTags, m_dimension);

        std::size_t totalElementHDNumber = 0;
        for(auto elmTags : dummyElementTags)
            totalElementHDNumber += elmTags.size();

        m_elementsHD.resize(totalElementHDNumber);

        std::size_t totalFaceNumber = 0;
        for(Entity& entity : m_entitiesHD)
        {
            entity.subTag = gmsh::model::addDiscreteEntity(m_dimension - 1);
            std::vector<std::size_t> nodesTagPerFace;
            gmsh::model::mesh::getElementEdgeNodes(entity.pElementProperty->type, nodesTagPerFace, entity.mainTag);
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
                    entity.pFaceProperty = &elmProp;
                    break;
                }
            }

            assert(entity.pFaceProperty != nullptr);

            gmsh::model::mesh::addElementsByType(entity.subTag, eleTypeSubEntity, {}, nodesTagPerFace);

            totalFaceNumber += nodesTagPerFace.size()/entity.pFaceProperty->numNodes;
        }

        m_faces.resize(totalFaceNumber);

        /** Fill the elements HD vector with th required infos. **/

        std::size_t elmCounter = 0;
        std::size_t faceCounter = 0;

        for(Entity& entity : m_entitiesHD)
        {
            std::vector<std::size_t> nodesTagPerFace;
            gmsh::model::mesh::getElementEdgeNodes(entity.pElementProperty->type, nodesTagPerFace, entity.mainTag);

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
            gmsh::model::mesh::getJacobians(entity.pFaceProperty->type, entity.pFaceProperty->intPointsCoord,
                                            jacobiansLD, determinantsLD, dummyCoordLD, entity.subTag);

            std::vector<double> baryCenters;
            gmsh::model::mesh::getBarycenters(entity.pElementProperty->type, entity.mainTag, false, true, baryCenters);

            unsigned int nNodesElm = entity.pElementProperty->numNodes;
            unsigned int nElements = elementTags.size();
            unsigned int nNodesFace = entity.pFaceProperty->numNodes;
            unsigned int nGPHD = entity.pElementProperty->intPointsWeigth.size();
            unsigned int nGPLD = entity.pFaceProperty->intPointsWeigth.size();
            unsigned int nFacesPerElm = determinantsLD.size()/(nGPLD*nElements);

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

                for(unsigned int i = 0 ; i < nFacesPerElm ; ++i)
                {
                    Face face = {};

                    face.nodesTag = std::vector<std::size_t>(nodesTagPerFace.begin() + nNodesFace*nFacesPerElm*e + nNodesFace*i,
                                                             nodesTagPerFace.begin() + nNodesFace*nFacesPerElm*e + nNodesFace*(i + 1));

                    for(std::size_t nodeTag : face.nodesTag)
                    {
                        std::vector<double> coord, dummyParametricCoord;
                        gmsh::model::mesh::getNode(nodeTag, coord, dummyParametricCoord);
                        face.nodesCoord.push_back(coord);
                    }

                    face.determinant = std::vector<double>(determinantsLD.begin() + nFacesPerElm*nGPLD*e + nGPLD*i,
                                                           determinantsLD.begin() + nFacesPerElm*nGPLD*e + nGPLD*(i + 1));

                    face.jacobian = std::vector<double>(jacobiansLD.begin() + nFacesPerElm*9*nGPLD*e + 9*nGPLD*i,
                                                        jacobiansLD.begin() + nFacesPerElm*9*nGPLD*e + 9*nGPLD*(i + 1));

                    computeFaceNormal(face, elementBarycenters);

                    m_faces[faceCounter] = std::move(face);
                    elm.pFaces.push_back(&m_faces.back());

                    faceCounter++;
                }

                m_elementsHD[elmCounter] = std::move(elm);
                entity.pElements.push_back(&m_elementsHD.back());

                Element& finalElm = m_elementsHD.back();
                for(Face* pFace : finalElm.pFaces)
                {
                    pFace->pParentElementHD = &finalElm;
                }

                elmCounter++;
            }
        }
    }

    void Mesh::computeFaceNormal(Face& face, const std::vector<double>& elementBarycenter)
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
        double nx = face.nodesCoord[1][1] - face.nodesCoord[0][1];
        double ny = face.nodesCoord[0][0] - face.nodesCoord[1][0];
        double norm = sqrt(ny*ny + nx*nx);

        // unfortunately, nodes per edge in nodes vector are not always in
        // the same order (clockwise vs anticlockwise) => we need to check
        // the orientation
        double vx = elementBarycenter[0] - (face.nodesCoord[1][0] + face.nodesCoord[0][0])/2;
        double vy = elementBarycenter[1] - (face.nodesCoord[1][1] + face.nodesCoord[0][1])/2;

        if(nx*vx + ny*vy > 0)
        {
            nx = -nx;
            ny = -ny;
        }

        // normalize the normal components
        normal.push_back(nx/norm);
        normal.push_back(ny/norm);

        face.normal = std::move(normal);
    }

    void Mesh::checkIfFaceIsBoundary()
    {
        assert(m_elementsLD.size() != 0);

        //TO DO: OpenMP This ! As well as use pointers to elementLD to remove elements found
        for(Face& face : m_faces)
        {
            for(Element& elementLD : m_elementsLD)
            {
                if(std::is_permutation(face.nodesTag.begin(), face.nodesTag.end(), elementLD.nodesTag.begin()))
                {
                    face.pElementLD = &elementLD;
                    break;
                }
            }
        }
    }

    void Mesh::associateFaces()
    {
        std::vector<Face*> pFaces;

        for(Face& face : m_faces)
        {
            if(face.pElementLD == nullptr)
                pFaces.push_back(&face);
        }

        assert(pFaces.size() % 2 == 0);

        while(pFaces.size() != 0)
        {
            std::vector<std::size_t> currentFaceNodeTag = pFaces[0]->nodesTag;
            std::size_t indexInFront = 0;

            #pragma omp parallel default(shared)
            {
                #pragma omp for schedule(static)
                for(std::size_t i = 1 ; i < pFaces.size() ; ++i)
                {
                    if(std::is_permutation(pFaces[i]->nodesTag.begin(), pFaces[i]->nodesTag.end(),
                                           currentFaceNodeTag.begin()))
                    {
                        #pragma omp critical
                        {
                            indexInFront = i;
                        }
                        #pragma omp cancel for
                    }

                    #pragma omp cancellation point for
                }
            }

            assert(indexInFront != 0);

            pFaces[0]->pFaceInFront = pFaces[indexInFront];
            pFaces[indexInFront]->pFaceInFront = pFaces[0];

            pFaces.erase(pFaces.begin() + indexInFront);
            pFaces.erase(pFaces.begin());
        }
    }
}
