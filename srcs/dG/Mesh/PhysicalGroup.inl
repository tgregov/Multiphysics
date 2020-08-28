#include "PhysicalGroup.hpp"

namespace dG
{
    inline int PhysicalGroup::getTag() const noexcept
    {
        return tag;
    }

    inline  std::string PhysicalGroup::getName() const noexcept
    {
        return name;
    }

    inline std::size_t PhysicalGroup::getEntityCount() const noexcept
    {
        return pEntities.size();
    }

    inline const Entity& PhysicalGroup::getEntity(std::size_t entity) const noexcept
    {
        return *pEntities[entity];
    }
}
