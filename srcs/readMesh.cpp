#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI
#include <gmsh.h>
#include "readMesh.hpp"


int posInVector(const std::vector<int>& vec, const std::pair<int, int>& couple)
{
    for(std::size_t i = 0; i < vec.size()/2; i++)
    {
        if(vec[2*i] == couple.first && vec[2*i+1] == couple.second)
            return i;
        if(vec[2*i] == couple.second && vec[2*i+1] == couple.first)
            return i;
    }

    return -1;
}


// TO DO: optimize this procedure (linked to the normal computation)
int diffVector(const std::vector<int>& vec, const std::pair<int, int>& couple)
{
    for(std::size_t i = 0; i < vec.size(); i++){
        if(vec[i] != couple.first && vec[i] != couple.second)
            return vec[i];
    }

    return -1;
}


bool readMesh(MeshParams& meshParams, const std::string& fileName, 
                const std::string& intScheme, const std::string& basisFuncType,
                const std::string& basisFuncGradType)
{
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
    int numComp;
    gmsh::model::mesh::getBasisFunctions(eleType2D, intScheme, basisFuncType,
                                         meshParams.intPoints, numComp, 
                                         meshParams.basisFunc);
    gmsh::model::mesh::getBasisFunctions(eleType2D, intScheme, basisFuncGradType,
                                         meshParams.intPoints, numComp, 
                                         meshParams.basisFuncGrad);


    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2);
    int c = entities[0].second; // c is the tag of the surface

    // Get the Jacobians information for the 2D triangular elements
    std::vector<double> jac, pts;
    gmsh::model::mesh::getJacobians(eleType2D, intScheme, meshParams.jacobian, 
                                    meshParams.determinant, pts, c);
    gmsh::finalize();

    meshParams.nGP = meshParams.intPoints.size()/4;
    meshParams.nSF = meshParams.basisFunc.size()/meshParams.nGP;
    meshParams.nE = meshParams.determinant.size()/meshParams.nGP;

    return true;
}
