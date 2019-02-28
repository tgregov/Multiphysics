#include <iostream>
#include <string>
#include <Eigen/Sparse>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
// #include <gmsh.h>
#include "readMesh.hpp"
#include "buildM.hpp"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 1;
    }

    MeshParams meshParams;

    if(!readMesh(meshParams, std::string(argv[1]), "Gauss1", "Lagrange"))
    {
        std::cerr << "[FAIL] The mesh was not read !" << std::endl;
        return -1;
    }
    else
    {
        std::cout << "The mesh was read successfully" << std::endl;
    }

    Eigen::SparseMatrix<double> M(meshParams.nE*meshParams.nSF,
                                    meshParams.nE*meshParams.nSF);

    buildM(meshParams, M);
    std::cout << std::endl;
    std::cout << M << std::endl;

    return 0;
}
