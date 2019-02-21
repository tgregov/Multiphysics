#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
// #include <gmsh.h>
#include "readMesh.hpp"
#include <Eigen/Dense>

int main(int argc, char **argv)
{
    Eigen::MatrixXd m2();
    Mesh* mesh = readMesh(argc, argv);

    if(mesh == nullptr){
        std::cout << "[FAIL] The mesh was not read !" << std::endl;
        return 1;
    } else{
        std::cout << "The mesh was read successfully" << std::endl;
    }
	
    return 0;
}
