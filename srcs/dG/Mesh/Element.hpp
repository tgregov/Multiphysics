#pragma once
#ifndef ELEMENT_HPP_INCLUDED
#define ELEMENT_HPP_INCLUDED

#include <vector>

#include "dgMesh_export.h"

namespace dG
{
    class Entity;
    class Face;
    class Mesh;

    class DG_MESH_API Element
    {
        public:
            Element()                                    = default;
            Element(const Element& element)              = default;
            Element& operator=(const Element& element)   = default;
            Element(Element&& entity)                    = default;
            Element& operator=(Element&& element)        = default;
            ~Element()                                   = default;

            inline std::size_t getTag() const noexcept;
            inline std::size_t getNodeTag(unsigned short nodeIndex) const noexcept;
            inline std::vector<double> getNodeCoord(unsigned short nodeIndex) const noexcept;
            inline double getDeterminant(unsigned short gaussPoint) const noexcept;
            inline double getJacobian(unsigned short gaussPoint, unsigned short i, unsigned short j) const noexcept;

            inline const Entity& getEntity() const noexcept;
            inline const Face& getFace(unsigned short faceIndex) const noexcept;

        private:
            std::size_t tag                                         = -1;
            std::vector<std::size_t> nodesTag                       = {};
            std::vector<std::vector<double>> nodesCoord             = {};

            std::vector<double> determinant                         = {};
            std::vector<double> jacobian                            = {};

            Entity* pEntity                                         = nullptr;
            std::vector<Face*> pFaces                               = {};

            friend class Mesh;
    };
}

#include "Element.inl"

#endif // ELEMENT_HPP_INCLUDED

