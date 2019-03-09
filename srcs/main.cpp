#include <iostream>
#include <string>
#include <Eigen/Sparse>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
// #include <gmsh.h>
#include "Mesh2D.hpp"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 1;
    }

    Mesh2D mesh;

    if(!readMesh2D(mesh, std::string(argv[1]), "Gauss1", "Lagrange", "GradLagrange"))
    {
        std::cerr<<"Something went wrong when reading mesh file: "<<argv[1]<<std::endl;
        return -1;
    }

    return 0;
}
