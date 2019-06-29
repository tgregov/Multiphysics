/**
 * \file displayMesh.cpp
 * \brief Implementation of the required function to display a Mesh2D struct.
 */

#include <iostream>
#include <map>
#include "Mesh.hpp"


// see .hpp for description
void displayMesh(const Mesh& mesh)
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
    // general information about the current entity
    std::cout 	<< "[Entity (" << 1 << ")]:\n"
                << "\t- Tag of the " << mesh.dim << "D entity: "
                << mesh.entityTagHD << "\n"
                << "\t- Tag of the " << mesh.dim-1 << "D entity: "
                << mesh.entityTagLD << "\n"
                << "\t- Number of " << mesh.dim << "D elements: "
                <<  mesh.elements.size() << "\n";

    // display the information about each element
    for(size_t j = 0 ; j < mesh.elements.size() ; ++j)
    {

        // display the information about the current element
        Element element = mesh.elements[j];
        std::cout 	<< "\t[Element (" << j << ")]:\n"
                    << "\t\t- Tag: " << element.elementTag << "\n"
                    << "\t\t- " << mesh.dim << "D type: "
                    << element.elementTypeHD << "\n"
                    << "\t\t- " << mesh.dim-1 << "D type: "
                    << element.elementTypeLD << "\n"
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
        for(size_t k = 0 ; k < element.determinantHD.size() ; ++k)
        {

            std::cout	<< "\t\t- [at GP (" << k << ")]:\n"
                        << "\t\t\t- det = " << element.determinantHD[k]
                        << std::endl;
            for(unsigned int l = 0 ; l < 9 ; ++l)
            {
                std::cout 	<< "\t\t\t- jac[" << l << "] = "
                            << element.jacobianHD[9*k + l] << std::endl;
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
            std::cout   << "\t\t\t- Normal: (" ;
            for(size_t y = 0 ;  y < edges[k].normal.size() ; ++y)
            {
                std::cout << element.edges[k].normal[y] << ", ";
            }
            std::cout   << ")\n"
                        << "\t\t\t- Det: ";
            for(size_t r = 0 ; r < edges[k].determinantLD.size() ; ++r)
            {
                std::cout << edges[k].determinantLD[r];
                if(r != edges[k].determinantLD.size() - 1)
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

                //Check that the normals are computed correctly
                for(unsigned int v = 0 ; v < edges[k].normal.size() ; ++v)
                {
                    if(edges[k].normal[v] != -mesh
                        .elements[edges[k].edgeInFront.first]
                        .edges[edges[k].edgeInFront.second].normal[v])
                    {
                        std::cerr 	<< "Bug in the normal of that edge !"
                                    << std::endl;
                        return;
                    }
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
	std::cout << std::endl;


	/*******************************************************************************
	 *			   DISPLAY INFORMATION ABOUT THE ELEMENT TYPES                     *
	 *******************************************************************************/
	// display the number of 2D element types
	std::cout 	<< "Number of element types: " << mesh.elementProperties.size()
				<< std::endl;

	// loop over the element types
	for(std::pair<int, ElementProperty> elmProperty : mesh.elementProperties)
	{
		std::cout 	<< "[Type (" << elmProperty.first << ")]:" << std::endl;

		ElementProperty property  = elmProperty.second;
		std::cout 	<< "\t- name: " 	<< property.name 		<< std::endl
					<< "\t- dim: " 		<< property.dim 		<< std::endl
					<< "\t- order: " 	<< property.order 		<< std::endl
					<< "\t- numNodes: " << property.numNodes 	<< std::endl;
	}
}
