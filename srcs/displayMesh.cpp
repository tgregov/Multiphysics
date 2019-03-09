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
					<< "	- Tag of the 2D entity: " << entitity.entityTag2D 
					<< std::endl
					<< "	- Tag of the 1D entity: " << entitity.entityTag1D
					<< std::endl;

		std::cout 	<< "	- Number of 2D elements: " 
					<<  entitity.elements.size() << std::endl;

		for(unsigned int j = 0 ; j < entitity.elements.size() ; ++j)
		{

			Element2D element = mesh.entities[i].elements[j];
			std::cout 	<< "	[Element (" << j << ")]:" << std::endl
						<< "		- Tag: " << element.elementTag << std::endl
			 			<< "		- 2D type: " << element.elementType2D 
						<< std::endl
						<< "		- 1D type: " << element.elementType1D 
						<< std::endl; 

			for(unsigned int k = 0 ; k < element.determinant2D.size() ; ++k)
			{
				std::cout	<< "		- [at GP (" << k << ")]:" << std::endl
							<< "			- det = " << element.determinant2D[k] 
							<< std::endl;
				for(unsigned int l = 0 ; l < 9 ; ++l)
				{
					std::cout 	<< "			- jac[" << l << "] = " 
								<< element.jacobian2D[9*k + l] << std::endl;
				}

			}

			std::vector<Edge> edges = element.edges;
			for(unsigned int k = 0 ; k < edges.size() ; ++k)
			{
				std::cout	<< "		- [Edge (" << k << ")]:" << std::endl
							<< "			- Tag A: " << edges[k].nodeTags.first 
							<< std::endl
							<< "			- Tag B: " << edges[k].nodeTags.second 
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

		std::cout 	<< "	- name: " << property.name << std::endl
					<< "	- dim: " << property.dim << std::endl	
					<< "	- order: " << property.order << std::endl	
					<< "	- numNodes: " << property.numNodes << std::endl;	
	}

	std::cout 	<< "Number of 1D element types: " << mesh.elementProperties1D.size() 
				<< std::endl;	

	for(std::pair<int, ElementProperty> elmProperty : mesh.elementProperties1D) 
	{
		std::cout 	<< "[Type (" << elmProperty.first << ")]:" << std::endl;
		ElementProperty property  = elmProperty.second; 

		std::cout 	<< "	- name: " << property.name << std::endl
					<< "	- dim: " << property.dim << std::endl	
					<< "	- order: " << property.order << std::endl	
					<< "	- numNodes: " << property.numNodes << std::endl;	
	}
}