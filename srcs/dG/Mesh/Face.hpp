#pragma once
#ifndef FACE_HPP_INCLUDED
#define FACE_HPP_INCLUDED

#include <vector>

#include "dgMesh_export.h"

namespace dG
{
    class Element;
    class Mesh;

    class DG_MESH_API Face
    {
        public:
            Face()                              = default;
            Face(const Face& face)              = default;
            Face& operator=(const Face& face)   = default;
            Face(Face&& face)                   = default;
            Face& operator=(Face&& face)        = default;
            ~Face()                             = default;

            inline std::size_t getNodeTag(unsigned short nodeIndex) const noexcept;
            inline std::vector<double> getNodeCoord(unsigned short nodeIndex) const noexcept;

            inline double getDeterminant(unsigned short gaussPoint) const noexcept;
            inline double getJacobian(unsigned short gaussPoint, unsigned short i, unsigned short j) const noexcept;

            inline std::vector<double> getNormal() const noexcept;

            inline const Element& getParentElement() const noexcept;
            inline bool isOnBOundary() const noexcept;
            inline const Element& getBoundaryElement() const noexcept;
            inline const Face& getFaceInFront() const noexcept;

        private:
            std::vector<std::size_t> nodesTag                       = {};
            std::vector<std::vector<double>> nodesCoord             = {};

            std::vector<double> determinant                         = {};
            std::vector<double> jacobian                            = {};

            std::vector<double> normal                              = {};

            Element* pParentElementHD                               = nullptr;
            Element* pElementLD                                     = nullptr;
            Face* pFaceInFront                                      = nullptr;

            friend class Mesh;
    };
}

#include "Face.inl"

#endif // FACE_HPP_INCLUDED
