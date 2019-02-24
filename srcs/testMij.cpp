#include <iostream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI
#include <gmsh.h>

bool testMij(const std::string& fileName)
{
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

    // -------------------------------------------------------------------------
    // 							GET BASIS FUNCTIONS PART 
    // -------------------------------------------------------------------------
    std::cout << "[GET BASIS FUNCTIONS]" << std::endl;
	gmsh::model::mesh::getElementTypes(eleTypes, 2);
    int eleType2D = eleTypes[0];
    std::cout << "Element type of dim = 2: " << eleType2D << std::endl;
    
    std::vector<double> intpts, bf;
    int numComp;
    gmsh::model::mesh::getBasisFunctions(eleType2D, "Gauss1", "Lagrange",
                                         intpts, numComp, bf);
    std::cout << "For Gauss 1, Lagrange:" << std::endl;
    std::cout << " - number of components: " << numComp << std::endl;
    std::cout << " - intpts: " << std::endl;
    std::cout << "    - length: " << intpts.size() << std::endl;
    for(size_t i = 0; i < intpts.size(); i++)
    {
    	std::cout << "    - intpts[" << i << "] = " << intpts[i] << std::endl;	
    }
   	std::cout << " - bf: " << std::endl;
    std::cout << "    - length: " << bf.size() << std::endl;
    for(size_t i = 0; i < bf.size(); i++)
    {
    	std::cout << "    - bf[" << i << "] = " << bf[i] << std::endl;	
    }

    // The way I understand it:
    // - intpts contains, for each GP points exactly 4 values: the (x, y, z) 
    //   coordinates + the associated weigth => length = 4*(number of GP)
    // - bf contains the evaluation of each GP, for each "component" of the 
    //	 basis function. For example in 2D, for Lagrange we have 3 "components", 
    //   given by the fact that there are 3 shape functions: 1 associated to
    //   corner. bf then contains, for each shape function its evaluation at all
    //	 the GP.
    // - numComp always return 1 => maybe a problem ? If I understand correctly
    //   what is meant by "component", it should not be like that

    // -------------------------------------------------------------------------
    // 							GET JACOBIAN PART 
    // -------------------------------------------------------------------------
    std::cout << "[GET JACOBIAN]" << std::endl;

    // Get the 2D entities => the surface on which the elements are defined
	std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2);
    int c = entities[0].second; // c is the tag of the surface 

    // Get the Jacobians information for the 2D triangular elements
    std::vector<double> jac, det, pts;
    gmsh::model::mesh::getJacobians(eleType2D, "Gauss1", jac, det, pts, c);
    std::cout << "For Gauss 1:" << std::endl;
    std::cout << " - jac:" << std::endl;
    std::cout << "    - length: " << jac.size() << std::endl;
    for(size_t i = 0; i < jac.size(); i++)
    {
    	std::cout << "    - jac[" << i << "] = " << jac[i] << std::endl;	
    }
    std::cout << " - det:" << std::endl;
    std::cout << "    - length: " << det.size() << std::endl;
    for(size_t i = 0; i < det.size(); i++)
    {
    	std::cout << "    - det[" << i << "] = " << det[i] << std::endl;	
    }
    std::cout << " - pts:" << std::endl;
    std::cout << "    - length: " << pts.size() << std::endl;
    for(size_t i = 0; i < pts.size(); i++)
    {
    	std::cout << "    - pts[" << i << "] = " << pts[i] << std::endl;	
    }

    // The way I understand it:
    // - jac contains for each element, for each GP, the 9 components of the 
    //   Jacobian matrix => in our case, we have for instance 1 GP, 4 elements
    //	 => we have 4*1*9 = 36 components for the vector jac
    // - det contains for each element, for each GP the determinant of the 
    //   Jacobian matrix => we have 4*1 = 4 components for the vector det
    // - pts contains for each element the (x, y, z) coordinates in the 
    //   physical frame of the GP => we have 4*3 = 12 components


    return true;
}



