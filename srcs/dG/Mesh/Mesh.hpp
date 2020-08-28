#pragma once
#ifndef MESH_HPP_INCLUDED
#define MESH_HPP_INCLUDED

#include <array>
#include <string>
#include <vector>

#include "PhysicalGroup.hpp"
#include "Entity.hpp"
#include "Element.hpp"
#include "Face.hpp"
#include "ElementProperty.hpp"

#include "dgMesh_export.h"

namespace dG
{
    class DG_MESH_API Mesh
    {
        public:
            Mesh()                              = delete;
            Mesh(std::string integrationType, std::string basisFuncType);
            Mesh(const Mesh& mesh)              = delete;
            Mesh& operator=(const Mesh& mesh)   = delete;
            Mesh(Mesh&& mesh)                   = delete;
            Mesh& operator=(Mesh&& mesh)        = delete;
            ~Mesh()                             = default;

            void displayToConsole() const noexcept;
            void loadFromFile(std::string fileName);

            inline unsigned int getEntityHDCount() const noexcept;
            inline const Entity& getEntityHD(unsigned int entityIndex) const noexcept;

        private:
            void loadPhysicalGroupsAndEntities();
            void loadElementsProperty();
            void loadElements();
            void checkIfFaceIsBoundary();
            void associateFaces();

            void computeFaceNormal(Face& face, const std::vector<double>& elementBarycenter);

            int m_dimension;

            std::string m_integrationType;
            std::string m_basisFuncType;

            std::vector<PhysicalGroup> m_physicalGroupHD;
            std::vector<PhysicalGroup> m_physicalGroupLD;

            std::vector<Entity> m_entitiesHD;
            std::vector<Entity> m_entitiesLD;

            std::vector<ElementProperty> m_elementsPropertyHD;
            std::vector<ElementProperty> m_elementsPropertyLD;

            std::vector<Element>  m_elementsHD;
            std::vector<Element>  m_elementsLD;
            std::vector<Face> m_faces;
    };
}

#include "Mesh.inl"

#endif // MESH_HPP_INCLUDED
