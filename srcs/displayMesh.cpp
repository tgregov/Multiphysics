/**
 * \file dispalyMesh.cpp
 * \brief Implementation of the required function to display a Mesh2D struct.
 */

#include <iostream>
#include <map>
#include "Mesh2D.hpp"


void displayMesh(const Mesh2D& mesh)
{
	// display the information about the entities
	std::cout 	<< "Number of entites in the mesh: " << mesh.entities.size()
				<< std::endl;

	for(unsigned int i = 0 ; i < mesh.entities.size() ; ++i)
	{
		Entity2D entitity =  mesh.entities[i];

		std::cout 	<< "[Entity (" << i << ")]:" << std::endl
					<< "\t- Tag of the 2D entity: " << entitity.entityTag2D 
					<< std::endl
					<< "\t- Tag of the 1D entity: " << entitity.entityTag1D
					<< std::endl;

		std::cout 	<< "\t- Number of 2D elements: "
					<<  entitity.elements.size() << std::endl;

		for(unsigned int j = 0 ; j < entitity.elements.size() ; ++j)
		{

			Element2D element = mesh.entities[i].elements[j];
			std::cout 	<< "\t[Element (" << j << ")]:" << std::endl
						<< "\t\t- Tag: " << element.elementTag << std::endl
			 			<< "\t\t- 2D type: " << element.elementType2D << std::endl
						<< "\t\t- 1D type: " << element.elementType1D << std::endl;

			for(unsigned int k = 0 ; k < element.determinant2D.size() ; ++k)
			{
				std::cout	<< "\t\t- [at GP (" << k << ")]:" << std::endl
							<< "\t\t\t- det = " << element.determinant2D[k] 
							<< std::endl;
				for(unsigned int l = 0 ; l < 9 ; ++l)
				{
					std::cout 	<< "\t\t\t- jac[" << l << "] = "
								<< element.jacobian2D[9*k + l] << std::endl;
				}

			}

			std::vector<Edge> edges = element.edges;
			for(unsigned int k = 0 ; k < edges.size() ; ++k)
			{
				std::cout	<< "\t\t- [Edge (" << k << ")]:" << std::endl
							<< "\t\t\t- Tag A: " << edges[k].nodeTags.first 
							<< std::endl
							<< "\t\t\t- Tag B: " << edges[k].nodeTags.second
							<< std::endl
							<< "\t\t\t- Normal: (" << element.edgesNormal[k].first
							<< ", " << element.edgesNormal[k].second << ")"
							<< std::endl
							<< "\t\t\t- Det: " << edges[k].determinant1D[0] 
							<< std::endl;
			}
		}
	}
	std::cout << std::endl;

	// display the information of the element types
	std::cout 	<< "Number of 2D element types: " << mesh.elementProperties2D.size()
				<< std::endl;

	for(std::pair<int, ElementProperty> elmProperty : mesh.elementProperties2D)
	{
		std::cout 	<< "[Type (" << elmProperty.first << ")]:" << std::endl;
		ElementProperty property  = elmProperty.second;

		std::cout 	<< "\t- name: " << property.name << std::endl
					<< "\t- dim: " << property.dim << std::endl
					<< "\t- order: " << property.order << std::endl
					<< "\t- numNodes: " << property.numNodes << std::endl;
	}

	std::cout 	<< "Number of 1D element types: " << mesh.elementProperties1D.size()
				<< std::endl;

	for(std::pair<int, ElementProperty> elmProperty : mesh.elementProperties1D)
	{
		std::cout 	<< "[Type (" << elmProperty.first << ")]:" << std::endl;
		ElementProperty property  = elmProperty.second;

		std::cout 	<< "\t- name: " << property.name << std::endl
					<< "\t- dim: " << property.dim << std::endl
					<< "\t- order: " << property.order << std::endl
					<< "\t- numNodes: " << property.numNodes << std::endl;
	}
}
