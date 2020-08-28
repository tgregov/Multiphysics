#include "Element.hpp"

namespace dG
{
    inline std::size_t Element::getTag() const noexcept
    {
        return tag;
    }

    inline std::size_t Element::getNodeTag(unsigned short nodeIndex) const noexcept
    {
        return nodesTag[nodeIndex];
    }

    inline std::vector<double> Element::getNodeCoord(unsigned short nodeIndex) const noexcept
    {
        return nodesCoord[nodeIndex];
    }

    inline double Element::getDeterminant(unsigned short gaussPoint) const noexcept
    {
        return determinant[gaussPoint];
    }

    inline double Element::getJacobian(unsigned short gaussPoint, unsigned short i, unsigned short j) const noexcept
    {
        return jacobian[9*gaussPoint + 3*i + j];
    }

    inline const Entity& Element::getEntity() const noexcept
    {
        return *pEntity;
    }

    inline const Face& Element::getFace(unsigned short faceIndex) const noexcept
    {
        return *pFaces[faceIndex];
    }
}
