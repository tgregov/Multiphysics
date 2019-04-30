#include <cmath>
#include <iostream>
#include <gmsh.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "computeL2Norm.hpp"
#include "computeError.hpp"


bool computeError(const Mesh& mesh, const std::string& meshName, const std::string& resultsName)
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

	double t = 0.03;
	int step = 300;

    //get the results computed by the solver
    gmsh::view::getModelData(tags[0], step, dataType, elementTags,
                data, t, numComponents);

    int totalSize = data.size()*data[0].size();

    //put the field in a vector
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
double error = computeL2Norm(mesh, t, u);
gmsh::finalize();

std::cout << "Error value:" << error << std::endl;
        

    return true;
}