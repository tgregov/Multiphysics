#pragma once
#ifndef ENTITY_HPP_INCLUDED
#define ENTITY_HPP_INCLUDED

#include <vector>

#include "dgMesh_export.h"

namespace dG
{
    class PhysicalGroup;
    struct ElementProperty;
    class Element;
    class Mesh;

    class DG_MESH_API Entity
    {
        public:
            Entity()                                  = default;
            Entity(const Entity& entity)              = default;
            Entity& operator=(const Entity& entity)   = default;
            Entity(Entity&& entity)                   = default;
            Entity& operator=(Entity&& entity)        = default;
            ~Entity()                                 = default;

            inline int getMainTag() const noexcept;
            inline int getSubTag() const noexcept;

            inline std::size_t getPhysicalGroupCount() const noexcept;
            inline const PhysicalGroup& getPhysicalGroup(std::size_t pg) const noexcept;

            inline const ElementProperty& getElementProperty() const noexcept;
            inline const ElementProperty& getFaceProperty() const noexcept;

            inline std::size_t getElementCount() const noexcept;
            inline const Element& getElement(std::size_t elmIndex) const noexcept;

        private:
            int mainTag                                  = -1;
            int subTag                                   = -1;
            std::vector<PhysicalGroup*> pPhysicalGroups  = {};

            ElementProperty* pElementProperty            = nullptr;
            ElementProperty* pFaceProperty               = nullptr;
            std::vector<Element*> pElements              = {};

            friend class Mesh;
    };
}

#include "Entity.inl"

#endif // ENTITY_HPP_INCLUDED

