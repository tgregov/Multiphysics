#pragma once
#ifndef PHYSICALGROUP_HPP_INCLUDED
#define PHYSICALGROUP_HPP_INCLUDED

#include <string>
#include <vector>

#include "dgMesh_export.h"

namespace dG
{
    class Entity;
    class Mesh;

    class DG_MESH_API PhysicalGroup
    {
        public:
            PhysicalGroup()                                         = default;
            PhysicalGroup(const PhysicalGroup& entity)              = default;
            PhysicalGroup& operator=(const PhysicalGroup& entity)   = default;
            PhysicalGroup(PhysicalGroup&& entity)                   = default;
            PhysicalGroup& operator=(PhysicalGroup&& entity)        = default;
            ~PhysicalGroup()                                        = default;

            inline int getTag() const noexcept;
            inline  std::string getName() const noexcept;

            inline std::size_t getEntityCount() const noexcept;
            inline const Entity& getEntity(std::size_t entity) const noexcept;

        private:
            int tag                         = -1;
            std::string name                = {};
            std::vector<Entity*> pEntities  = {};

            friend class Mesh;
    };
}

#include "PhysicalGroup.inl"

#endif // PHYSICALGROUP_HPP_INCLUDED
