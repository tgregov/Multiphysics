#include "Mesh.hpp"

namespace dG
{
    inline unsigned int Mesh::getEntityHDCount() const noexcept
    {
        return m_entitiesHD.size();
    }

    inline const Entity& Mesh::getEntityHD(unsigned int entityIndex) const noexcept
    {
        return m_entitiesHD[entityIndex];
    }
}
