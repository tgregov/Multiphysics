#include <iostream>
#include <string>
#include <Eigen/Sparse>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
// #include <gmsh.h>
#include "Mesh2D.hpp"
#include "buildM.hpp"
#include "buildS.hpp"
#include "displayMesh.hpp"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 1;
    }

    Mesh2D mesh;

    if(!readMesh2D(mesh, std::string(argv[1]), "Gauss1", "Lagrange"))
    {
        std::cerr   << "Something went wrong when reading mesh file: "
                    << argv[1] << std::endl;
        return -1;
    }

    displayMesh(mesh);

    unsigned long numNodes = getNumNodes(mesh);
    Eigen::SparseMatrix<double> M(numNodes, numNodes);
    Eigen::SparseMatrix<double> Sx(numNodes, numNodes);
    Eigen::SparseMatrix<double> Sy(numNodes, numNodes);
    buildM(mesh, M);
    buildS(mesh, Sx, Sy);
    std::cout << "Matrix [M]:\n" << M << std::endl;
    std::cout << "Matrix [Sx]:\n" << Sx << std::endl;
    std::cout << "Matrix [Sy]:\n" << Sy << std::endl;


    return 0;
}
