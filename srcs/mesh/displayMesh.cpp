/**
 * \file displayMesh.cpp
 * \brief Implementation of the required function to display a Mesh2D struct.
 */

#include <iostream>
#include <map>
#include "Mesh2D.hpp"


void displayMesh(const Mesh2D& mesh)
{

	// display the number of entities
	std::cout 	<< "Number of entites in the mesh: " << mesh.entities.size()
				<< std::endl;

	// display the information of each entity
	for(unsigned int i = 0 ; i < mesh.entities.size() ; ++i)
	{

		Entity2D entitity =  mesh.entities[i];

		// general information about the current entity
		std::cout 	<< "[Entity (" << i << ")]:" << std::endl
					<< "\t- Tag of the 2D entity: " << entitity.entityTag2D
					<< std::endl
					<< "\t- Tag of the 1D entity: " << entitity.entityTag1D
					<< std::endl;

		std::cout 	<< "\t- Number of 2D elements: "
					<<  entitity.elements.size() << std::endl;

		// display the information about each element
		for(unsigned int j = 0 ; j < entitity.elements.size() ; ++j)
		{

			// display the information about the current element
			Element2D element = mesh.entities[i].elements[j];
			std::cout 	<< "\t[Element (" << j << ")]:" << std::endl
						<< "\t\t- Tag: " << element.elementTag << std::endl
			 			<< "\t\t- 2D type: " << element.elementType2D << std::endl
						<< "\t\t- 1D type: " << element.elementType1D << std::endl
						<< "\t\t- Offset in u: " << element.offsetInU << std::endl
						<< "\t\t- Nodes tag: ";
            for(unsigned int p = 0 ; p < element.nodeTags.size() ; ++p)
                std::cout<<element.nodeTags[p]<<", ";

            std::cout   << std::endl;

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

			// display the information about the edges
			std::vector<Edge> edges = element.edges;
			for(unsigned int k = 0 ; k < edges.size() ; ++k)
			{
				std::cout	<< "\t\t- [Edge (" << k << ")]:" << std::endl
							<< "\t\t\t- Tag A: " << edges[k].nodeTags[0]
							<< std::endl
							<< "\t\t\t- Tag B: " << edges[k].nodeTags[1]
							<< std::endl
							<< "\t\t\t- Normal: (" << element.edges[k].normal.first
							<< ", " << element.edges[k].normal.second << ")"
							<< std::endl
							<< "\t\t\t- Det: ";
				for(unsigned int r = 0 ; r < edges[k].determinant1D.size() ; ++r)
               	{
					std::cout << edges[k].determinant1D[r];
					if(r != edges[k].determinant1D.size() - 1)
					{
						std::cout << ", ";
					}
					else
					{
						std::cout << std::endl;
					}
				}

				for(unsigned int i = 0 ; i < edges[k].offsetInElm.size() ; ++i)
				{
				std::cout 	<< "\t\t\t- OffsetInElm of node (" << i << "): "
							<< edges[k].offsetInElm[i]
							<< std::endl;
				}

                if(edges[k].edgeInFront.first != -1)
                {
                    std::cout << "\t\t\t- Edge in front: " << "element "
                                << edges[k].edgeInFront.first<<", "
                                << "edge " << edges[k].edgeInFront.second << ", "
                                << std::endl; //[TO DO]: put again inverted
                    if((edges[k].normal.first != -entitity
                    	.elements[edges[k].edgeInFront.first]
                    	.edges[edges[k].edgeInFront.second].normal.first)
                       || (edges[k].normal.second != -entitity
                       	.elements[edges[k].edgeInFront.first].
                       	edges[edges[k].edgeInFront.second].normal.second))
                    {
                        std::cerr << "Bug in the normal of that edge !" << std::endl;
                        return;
                    }
                }
                else
                {
                    if(edges[k].bcName.size() == 0)
                    {
                        std::cerr 	<< "BC node does not have a BC name !"
                        			<< std::endl;
                        return;
                    }
                    else
                        std::cout << "\t\t\t- BC: " << edges[k].bcName << std::endl;
                }
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
