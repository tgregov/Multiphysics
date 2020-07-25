#include <iostream>
#include <string>

#include "dgMesh/dgMesh.hpp"

int main(int argc, char **argv)
{
    dgMesh mesh("Gauss1", "Lagrange");
    mesh.loadFromFile(std::string(argv[1]));
    mesh.displayToConsole();

    return 0;
}
