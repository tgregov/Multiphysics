/**
 * \file displayMesh.cpp
 * \brief Implementation of the required function to display a Mesh2D struct.
 */

#include <iostream>
#include <map>
#include "Mesh2D.hpp"


// see .hpp for description
void displayMesh(const Mesh2D& mesh)
{

    std::cout   << "================================================================"
                << std::endl
                << "                 DISPLAYING THE MESH INFORMATION                "
                << std::endl
                << "================================================================"
                << std::endl;

	/*******************************************************************************
	 *					DISPLAY INFORMATION ABOUT THE ENTITIES                     *
	 *******************************************************************************/
	// display the number of entities
	std::cout 	<< "Number of entites in the mesh: " << mesh.entities.size()
				<< std::endl;

	// display the information of each entity
	for(size_t i = 0 ; i < mesh.entities.size() ; ++i)
	{

		// current entity
		Entity2D entitity =  mesh.entities[i];

		// general information about the current entity
		std::cout 	<< "[Entity (" << i << ")]:\n" 	
					<< "\t- Tag of the 2D entity: " 
					<< entitity.entityTag2D << "\n"
					<< "\t- Tag of the 1D entity: " 
					<< entitity.entityTag1D << "\n"
					<< "\t- Number of 2D elements: "
					<<  entitity.elements.size() << "\n";

		// display the information about each element
		for(size_t j = 0 ; j < entitity.elements.size() ; ++j)
		{

			// display the information about the current element
			Element2D element = mesh.entities[i].elements[j];
			std::cout 	<< "\t[Element (" << j << ")]:\n"
						<< "\t\t- Tag: " << element.elementTag << "\n"
			 			<< "\t\t- 2D type: " << element.elementType2D << "\n"
						<< "\t\t- 1D type: " << element.elementType1D << "\n"
						<< "\t\t- Offset in u: " << element.offsetInU << "\n"
						<< "\t\t- Nodes tags: " << std::endl;

			// node tags
            for(size_t p = 0 ; p < element.nodeTags.size() ; ++p)
            {
                std::cout	<<	"\t\t\t Node (" << p << "): "
                			<<  element.nodeTags[p] << std::endl;
            }
            std::cout 	<< std::endl;

            // determinant of the change of variable
			for(size_t k = 0 ; k < element.determinant2D.size() ; ++k)
			{

				std::cout	<< "\t\t- [at GP (" << k << ")]:\n"
							<< "\t\t\t- det = " << element.determinant2D[k]
							<< std::endl;
				for(unsigned int l = 0 ; l < 9 ; ++l)
				{
					std::cout 	<< "\t\t\t- jac[" << l << "] = "
								<< element.jacobian2D[9*k + l] << std::endl;
				}
			}

			// edges
			std::vector<Edge> edges = element.edges;
			for(size_t k = 0 ; k < edges.size() ; ++k)
			{

				std::cout	<< "\t\t- [Edge (" << k << ")]:\n";

				// node tags
				for(size_t b = 0 ; b < edges[k].nodeTags.size() ; ++b)
                {
                        std::cout 	<< "\t\t\t- Tag (" << b << "): " 
                        			<< edges[k].nodeTags[b] << std::endl;
                }

                // normal and determinant
                std::cout   << "\t\t\t- Normal: (" << element.edges[k].normal.first
							<< ", " << element.edges[k].normal.second << ")\n"
							<< "\t\t\t- Det: ";
				for(size_t r = 0 ; r < edges[k].determinant1D.size() ; ++r)
               	{
					std::cout << edges[k].determinant1D[r];
					if(r != edges[k].determinant1D.size() - 1)
					{
						std::cout 	<< ", ";
					}
					else
					{
						std::cout 	<< std::endl;
					}
				}

				// offset in the element
				for(size_t i = 0 ; i < edges[k].offsetInElm.size() ; ++i)
				{
				std::cout 	<< "\t\t\t- OffsetInElm of node (" << i << "): "
							<< edges[k].offsetInElm[i]
							<< std::endl;
				}

				// neighbour or BC
                if(edges[k].edgeInFront.first != -1)
                {
                    std::cout 	<< "\t\t\t- Edge in front: element "
                                << edges[k].edgeInFront.first << ", "
                                << "edge " << edges[k].edgeInFront.second
                                << std::endl;

                   	// checl that the normals are computed correctly
                    if((edges[k].normal.first != -entitity
                    	.elements[edges[k].edgeInFront.first]
                    	.edges[edges[k].edgeInFront.second].normal.first)
                       || (edges[k].normal.second != -entitity
                       	.elements[edges[k].edgeInFront.first].
                       	edges[edges[k].edgeInFront.second].normal.second))
                    {
                        std::cerr 	<< "Bug in the normal of that edge !" 
                        			<< std::endl;
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
                        std::cout 	<< "\t\t\t- BC: " << edges[k].bcName 
                    				<< std::endl;
                }
			}
		}
	}
	std::cout << std::endl;


	/*******************************************************************************
	 *			   DISPLAY INFORMATION ABOUT THE ELEMENT TYPES                     *
	 *******************************************************************************/
	// display the number of 2D element types
	std::cout 	<< "Number of 2D element types: " << mesh.elementProperties2D.size()
				<< std::endl;

	// loop over the 2D element types
	for(std::pair<int, ElementProperty> elmProperty : mesh.elementProperties2D)
	{
		std::cout 	<< "[Type (" << elmProperty.first << ")]:" << std::endl;

		ElementProperty property  = elmProperty.second;
		std::cout 	<< "\t- name: " 	<< property.name 		<< std::endl
					<< "\t- dim: " 		<< property.dim 		<< std::endl
					<< "\t- order: " 	<< property.order 		<< std::endl
					<< "\t- numNodes: " << property.numNodes 	<< std::endl;
	}

	// display the number of 1D element types
	std::cout 	<< "Number of 1D element types: " << mesh.elementProperties1D.size()
				<< std::endl;

	// loop over the 1D element types
	for(std::pair<int, ElementProperty> elmProperty : mesh.elementProperties1D)
	{
		std::cout 	<< "[Type (" << elmProperty.first << ")]:" << std::endl;

		ElementProperty property  = elmProperty.second;
		std::cout 	<< "\t- name: " 	<< property.name 		<< std::endl
					<< "\t- dim: " 		<< property.dim 		<< std::endl
					<< "\t- order: " 	<< property.order 		<< std::endl
					<< "\t- numNodes: " << property.numNodes 	<< std::endl;
	}
}
