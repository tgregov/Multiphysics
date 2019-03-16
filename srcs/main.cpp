#include <iostream>
#include <string>
#include <Eigen/Sparse>
#include "Mesh2D.hpp"
#include "displayMesh.hpp"
#include "timeInteg.hpp"

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
    if(!timeInteg(mesh, "RK1", 1, 3, "strong", std::string(argv[1])))
    {
        std::cerr   << "Something went wrong when time integrating" << std::endl;
        return -1;        
    }

    return 0;
}
