#include "Entity.hpp"

namespace dG
{
    inline int Entity::getMainTag() const noexcept
    {
        return mainTag;
    }

    inline int Entity::getSubTag() const noexcept
    {
        return subTag;
    }

    inline std::size_t Entity::getPhysicalGroupCount() const noexcept
    {
        return pPhysicalGroups.size();
    }

    inline const PhysicalGroup& Entity::getPhysicalGroup(std::size_t pg) const noexcept
    {
        return *pPhysicalGroups[pg];
    }

    inline const ElementProperty& Entity::getElementProperty() const noexcept
    {
        return *pElementProperty;
    }

    inline const ElementProperty& Entity::getFaceProperty() const noexcept
    {
        return *pFaceProperty;
    }

    inline std::size_t Entity::getElementCount() const noexcept
    {
        return pElements.size();
    }

    inline const Element& Entity::getElement(std::size_t elmIndex) const noexcept
    {
        return *pElements[elmIndex];
    }
}
