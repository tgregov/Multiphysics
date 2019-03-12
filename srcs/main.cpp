#include <iostream>
#include <string>
#include <Eigen/Sparse>
#include "Mesh2D.hpp"
#include "buildM.hpp"
#include "buildS.hpp"
#include "displayMesh.hpp"
#include "buildFlux.hpp"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 1;
    }

    Mesh2D mesh;

    if(!readMesh2D(mesh, std::string(argv[1]), "Gauss3", "Lagrange"))
    {
        std::cerr   << "Something went wrong when reading mesh file: "
                    << argv[1] << std::endl;
        return -1;
    }

    displayMesh(mesh);

    unsigned long numNodes = getNumNodes(mesh);
    /*Eigen::SparseMatrix<double> M(numNodes, numNodes);
    Eigen::SparseMatrix<double> Sx(numNodes, numNodes);
    Eigen::SparseMatrix<double> Sy(numNodes, numNodes);
    buildM(mesh, M);
    buildS(mesh, Sx, Sy);
    std::cout << "Matrix [M]:\n" << M << std::endl;
    std::cout << "Matrix [Sx]:\n" << Sx << std::endl;
    std::cout << "Matrix [Sy]:\n" << Sy << std::endl;*/
    Eigen::VectorXd I(numNodes); I.setZero();
    Eigen::VectorXd u(numNodes); u.setZero();
    buildFlux(mesh, I, u, "weak", numNodes);
    std::cout << "Vector I: \n" << I << std::endl;


    return 0;
}
