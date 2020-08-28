#include "Mesh.hpp"

namespace dG
{
    inline std::size_t Mesh::getEntityHDCount() const noexcept
    {
        return m_entitiesHD.size();
    }

    inline const Entity& Mesh::getEntityHD(unsigned int entityIndex) const noexcept
    {
        return m_entitiesHD[entityIndex];
    }
}
