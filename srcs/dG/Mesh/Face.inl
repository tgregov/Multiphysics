#include "Face.hpp"

namespace dG
{
    inline std::size_t Face::getNodeTag(unsigned short nodeIndex) const noexcept
    {
        return nodesTag[nodeIndex];
    }

    inline std::vector<double> Face::getNodeCoord(unsigned short nodeIndex) const noexcept
    {
        return nodesCoord[nodeIndex];
    }

    inline double Face::getDeterminant(unsigned short gaussPoint) const noexcept
    {
        return determinant[gaussPoint];
    }

    inline double Face::getJacobian(unsigned short gaussPoint, unsigned short i, unsigned short j) const noexcept
    {
        return jacobian[9*gaussPoint + 3*i + j];
    }

    inline std::vector<double> Face::getNormal() const noexcept
    {
        return normal;
    }

    inline const Element& Face::getParentElement() const noexcept
    {
        return *pParentElementHD;
    }

    inline bool Face::isOnBOundary() const noexcept
    {
        return (pElementLD != nullptr);
    }

    inline const Element& Face::getBoundaryElement() const noexcept
    {
        return *pElementLD;
    }

    inline const Face& Face::getFaceInFront() const noexcept
    {
        return *pFaceInFront;
    }
}
