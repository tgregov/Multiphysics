#include <cmath>
#include <iostream>
#include <gmsh.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "computeNorm.hpp"
#include "computeNormBis.hpp"
#include "computeError.hpp"


bool computeError(const Mesh& mesh, const SolverParams& solverParams, const std::string& meshName,
                     const std::string& resultsName, double& errorL2, double& errorLinf)
{
	gmsh::initialize();
	gmsh::option::setNumber("General.Terminal", 1);
	gmsh::open(resultsName);
	std::vector<int> tags;
	gmsh::view::getTags(tags);

	std::string dataType;
	int numComponents;
	std::vector<int> elementTags;
	std::vector<unsigned int> elementNumNodes = mesh.nodeData.elementNumNodes;
    std::vector<std::vector<double>> data(elementNumNodes.size());

	//double t1 = solverParams.simTime;
    double t1 = 0;
    double t2;
	int step = t1/solverParams.timeStep;


    //get the results computed by the solver, at timestep "step"
    gmsh::view::getModelData(tags[0], step, dataType, elementTags,
                data, t2, numComponents);

    std::cout << "t1,t2,step:" << std::endl << t1 << std::endl << t2 << std::endl << step << std::endl;
    // if (t1 != t2)
    // {
    //     std::cerr << ("ERROR: timestep mismatching") << std::endl;
    //     return false;
    // }
    
    int totalSize = data.size()*data[0].size();

    //put the results in a vector
    Eigen::VectorXd u(totalSize);
    int k=0;
    for (int i = 0; i < data.size(); ++i)
    {
    	for (int j = 0; j < data[i].size(); ++j)
    	{
    		u[k] = data[i][j];
    		++k;
    	}
    }


gmsh::finalize();

gmsh::initialize();
gmsh::option::setNumber("General.Terminal", 1);
gmsh::open(meshName);
computeNormBis(mesh, solverParams, t2, u, errorL2, errorLinf);
gmsh::finalize();

std::cout << "Error value:" << errorL2  << "\t" << errorLinf << std::endl;
        

    return true;
}