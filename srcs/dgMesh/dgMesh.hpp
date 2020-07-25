#pragma once
#ifndef DGMESH_HPP_INCLUDED
#define DGMESH_HPP_INCLUDED

#include <array>
#include <string>
#include <vector>

#include "dgMesh_export.h"

class DG_MESH_API dgMesh
{
    public:
        dgMesh()                                = delete;
        dgMesh(std::string integrationType, std::string basisFuncType);
        dgMesh(const dgMesh& mesh)              = delete;
        dgMesh& operator=(const dgMesh& mesh)   = delete;
        dgMesh(dgMesh&& mesh)                   = delete;
        dgMesh& operator=(dgMesh&& mesh)        = delete;
        ~dgMesh()                               = default;

        void displayToConsole() const noexcept;
        void loadFromFile(std::string fileName);

    private:
        void loadPhysicalGroupsAndEntities();
        void loadElementsProperty();
        void loadElements();
        void loadFaces();

        int m_dimension;

        std::string m_integrationType;
        std::string m_basisFuncType;

        struct PhysicalGroup;
        struct Entity;

        struct ElementProperty
        {
            int type                            = -1;
            std::string name                    = {};
            int order                           = -1;
            int numNodes                        = -1;
            std::vector<double> localNodeCoord  = {};
            int numPrimaryNodes                 = -1;

            std::vector<double> intPointsCoord  = {};
            std::vector<double> intPointsWeigth = {};

            int basisFuncNumComp                = -1;
            std::vector<double> basisFunctions  = {};

            int basisFuncNumCompGrad                = -1;
            std::vector<double> basisFunctionsGrad  = {};
        };

        struct Element
        {
            int tag                                                 = -1;
            std::vector<int> nodesTag                               = {};
            std::vector<std::vector<double>> nodesCoord             = {};

            std::vector<double> determinant                         = {};
            std::vector<std::vector<std::vector<double>>> jacobian  = {};

            std::vector<Element*> pElements                         = {};
            std::vector<Entity*> pEntity                            = {};
        };

        struct Entity
        {
            int tag                                      = -1;
            std::vector<PhysicalGroup*> pPhysicalGroups  = {};

            ElementProperty* pElementProperty            = nullptr;
            std::vector<Element*> pElements;
        };

        struct PhysicalGroup
        {
            int tag                         = -1;
            std::string name                = {};
            std::vector<Entity*> pEntities  = {};
        };

        std::vector<PhysicalGroup> m_physicalGroupHD;
        std::vector<PhysicalGroup> m_physicalGroupLD;

        std::vector<Entity> m_entitiesHD;
        std::vector<Entity> m_entitiesLD;

        std::vector<ElementProperty> m_elementsPropertyHD;
        std::vector<ElementProperty> m_elementsPropertyLD;

        std::vector<Element> m_elementsHD;
        std::vector<Element> m_elementsLD;
};

#include "dgMesh.inl"

#endif // DGMESH_HPP_INCLUDED
