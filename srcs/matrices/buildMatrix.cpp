#include <iostream>
#include "buildM.hpp"
#include "buildS.hpp"
#include "buildMatrix.hpp"

// see .hpp file for description
void buildMatrix(const Mesh& mesh, Matrix& matrix,
                 const DomainDiv& domainDiv, unsigned int rank)
{
    // redimension the matrix sizes
    matrix.invM.resize(domainDiv.node[rank], domainDiv.node[rank]);
    matrix.Sx.resize(domainDiv.node[rank], domainDiv.node[rank]);
    matrix.Sy.resize(domainDiv.node[rank], domainDiv.node[rank]);

    // build the invM matrix
    if(rank == 0)
        std::cout   << "Building the invM matrix...";
    buildM(mesh, matrix.invM, domainDiv, rank);
    // std::cout << "invM:\n" << matrix.invM;
    if(rank == 0)
    {
        std::cout   << "\rBuilding the invM matrix...       Done"       << std::flush
                    << std::endl;
    }


    // build the Sx and Sy matrices
    if(rank == 0)
        std::cout   << "Building the Sx and Sy matrices...";
    buildS(mesh, matrix.Sx, matrix.Sy, domainDiv, rank);
    // std::cout << "Sx:\n" << matrix.Sx;
    // std::cout << "Sy:\n" << matrix.Sy;
    if(rank == 0)
    {
        std::cout   << "\rBuilding the Sx and Sy matrices...    Done"   << std::flush
                    << std::endl;
    }
}
