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

    Mesh mesh;

    if(!readMesh(mesh, std::string(argv[1]))){
        std::cerr << "[FAIL] The mesh was not read !" << std::endl;
        return 1;
    } else{
        std::cout << "The mesh was read successfully" << std::endl;
    }

    /*
     * TEMPORARY BEGIN
     */

    // check that a .msh file was introduced
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open(fileName);

    // get the properties of 2D elements
    std::vector<int> eleTypes;
    gmsh::model::mesh::getElementTypes(eleTypes, 2);
    if (eleTypes.size() != 1)
    {
        // TO DO: handle hybrid meshes
        gmsh::logger::write("Hybrid meshes not handled in this example!",
                            "error");

        gmsh::finalize();
        return false;
    }
    int eleType2D = eleTypes[0]; // e.g. T3 elements

    // get basis functions
    std::vector<double> intpts, bf;
    int numComp;
    gmsh::model::mesh::getBasisFunctions(eleType2D, intScheme, basisFunc,
                                         intpts, numComp, bf);

    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2);
    int c = entities[0].second; // c is the tag of the surface 

    // Get the Jacobians information for the 2D triangular elements
    std::vector<double> jac, det, pts;
    gmsh::model::mesh::getJacobians(eleType2D, intScheme, jac, det, pts, c);
    gmsh::finalize();
 
    /*
     * TEMPORARY END
     */
    Eigen::SparseMatrix<double> M(nE*nSF, nE*nSF);


    return 0;
}