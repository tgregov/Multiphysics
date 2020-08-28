#include <iostream>
#include <stdexcept>
#include <string>

#include "dG/Mesh/Mesh.hpp"

int main(int argc, char **argv)
{
    dG::Mesh mesh("Gauss1", "Lagrange");

    try
    {
        mesh.loadFromFile(std::string(argv[1]));
        mesh.displayToConsole();
    }
    catch(const std::exception& exception)
    {
        std::cerr << "Something went wrong: " << exception.what() << std::endl;
        return -1;
    }
    catch(int e)
    {
        std::cerr << "GMSH certainly thrown something (only one to use int exceptions): "
                  << e << std::endl;
        return -2;
    }
    catch(...)
    {
        std::cerr << "Unexpected exception: bip bip boup ?" << std::endl;
        return -3;
    }

    return 0;
}
