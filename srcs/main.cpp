#include <iostream>
#include <string>
#include <Eigen/Dense>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
// #include <gmsh.h>
#include "readMesh.hpp"
#include "testMij.hpp"

  
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " file.msh [options]" << std::endl;
        return 1;
    }

    Eigen::MatrixXd m2(10, 10);
    Mesh mesh;

    if(!readMesh(mesh, std::string(argv[1]))){
        std::cerr << "[FAIL] The mesh was not read !" << std::endl;
        return 1;
    } else{
        std::cout << "The mesh was read successfully" << std::endl;
    }

    // TEMPORARY
    testMij(std::string(argv[1]));

    return 0;
}