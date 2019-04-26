#include <iostream>
#include "buildM.hpp"
#include "buildS.hpp"
#include "buildMatrix.hpp"

// see .hpp file for description
void buildMatrix(const Mesh& mesh, Matrix& matrix)
{
    // redimension the matrix sizes
    matrix.invM.resize(mesh.nodeData.numNodes, mesh.nodeData.numNodes);
    matrix.Sx.resize(mesh.nodeData.numNodes, mesh.nodeData.numNodes);
    matrix.Sy.resize(mesh.nodeData.numNodes, mesh.nodeData.numNodes);

    // build the invM matrix
    std::cout   << "Building the invM matrix...";
    buildM(mesh, matrix.invM);
    // std::cout << "invM:\n" << matrix.invM;
    std::cout   << "\rBuilding the invM matrix...       Done"       << std::flush
                << std::endl;

    // build the Sx and Sy matrices
    std::cout   << "Building the Sx and Sy matrices...";
    buildS(mesh, matrix.Sx, matrix.Sy);
    // std::cout << "Sx:\n" << matrix.Sx;
    // std::cout << "Sy:\n" << matrix.Sy;
    std::cout   << "\rBuilding the Sx and Sy matrices...    Done"   << std::flush
                << std::endl;

}
