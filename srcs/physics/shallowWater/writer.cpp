#include <gmsh.h>
#include "writer.hpp"

void writeShallow(std::vector<std::vector<double>>& uDisplay,
                  const std::vector<unsigned int>& elementNumNodes,
                  const std::vector<std::size_t>& elementTags,
                  const std::string& modelName, unsigned int nbreStep, double t,
                  const Field& field, const std::vector<double>& fluxCoeffs,
                  const std::vector<bool>& whatToWrite, std::vector<int>& viewTags)
{
    if(nbreStep == 0)
    {
        if(whatToWrite[0] == true)
            viewTags[0] = gmsh::view::add("H");

        if(whatToWrite[1] == true)
            viewTags[1] = gmsh::view::add("u");

        if(whatToWrite[2] == true)
            viewTags[2] = gmsh::view::add("v");

        if(whatToWrite[3] == true)
            viewTags[3] = gmsh::view::add("Specific KE");

        if(whatToWrite[4] == true)
            viewTags[4] = gmsh::view::add("Velocity Field");
    }

    if(whatToWrite[0] == true)
    {
        unsigned int offset = 0;
        for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
        {
            std::vector<double> temp(elementNumNodes[count]);
            for (unsigned int countLocal = 0; countLocal < elementNumNodes[count];
                ++countLocal)
            {
                temp[countLocal] = field.u[0][countLocal+offset];
            }
            offset += elementNumNodes[count];
            uDisplay[count] = std::move(temp);
        }

        gmsh::view::addModelData(viewTags[0], nbreStep, modelName,
                                 "ElementNodeData", elementTags, uDisplay, t, 1);
    }

    if(whatToWrite[1] == true)
    {
        unsigned int offset = 0;
        for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
        {
            std::vector<double> temp(elementNumNodes[count]);
            for (unsigned int countLocal = 0; countLocal < elementNumNodes[count];
                ++countLocal)
            {
                temp[countLocal] = field.u[1][countLocal+offset]
                                    /field.u[0][countLocal+offset];
            }
            offset += elementNumNodes[count];
            uDisplay[count] = std::move(temp);
        }

        gmsh::view::addModelData(viewTags[1], nbreStep, modelName,
                                 "ElementNodeData", elementTags, uDisplay, t, 1);
    }

    if(whatToWrite[2] == true)
    {
        unsigned int offset = 0;
        for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
        {
            std::vector<double> temp(elementNumNodes[count]);
            for (unsigned int countLocal = 0; countLocal < elementNumNodes[count];
                ++countLocal)
            {
                temp[countLocal] = field.u[2][countLocal+offset]
                                    /field.u[0][countLocal+offset];
            }
            offset += elementNumNodes[count];
            uDisplay[count] = std::move(temp);
        }

        gmsh::view::addModelData(viewTags[2], nbreStep, modelName,
                                 "ElementNodeData", elementTags, uDisplay, t, 1);
    }

    if(whatToWrite[3] == true)
    {
        unsigned int offset = 0;
        for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
        {
            std::vector<double> temp(elementNumNodes[count]);
            for (unsigned int countLocal = 0; countLocal < elementNumNodes[count];
                ++countLocal)
            {
                temp[countLocal] =
                0.5*(field.u[1][countLocal+offset]*field.u[1][countLocal+offset]
                     +field.u[2][countLocal+offset]*field.u[2][countLocal+offset])
                     /(field.u[0][countLocal+offset]*field.u[0][countLocal+offset]);
            }
            offset += elementNumNodes[count];
            uDisplay[count] = std::move(temp);
        }

        gmsh::view::addModelData(viewTags[3], nbreStep, modelName,
                                 "ElementNodeData", elementTags, uDisplay, t, 1);
    }

    if(whatToWrite[4] == true)
    {
        unsigned int offset = 0;
        for(size_t count = 0 ; count < elementNumNodes.size() ; ++count)
        {
            std::vector<double> temp(3*elementNumNodes[count]);
            for (unsigned int countLocal = 0; countLocal < elementNumNodes[count];
                ++countLocal)
            {
                temp[3*countLocal] = field.u[1][countLocal+offset]
                                     /field.u[0][countLocal+offset];

                temp[3*countLocal+1] = field.u[2][countLocal+offset]
                                       /field.u[0][countLocal+offset];

                temp[3*countLocal+2] = 0;
            }
            offset += elementNumNodes[count];
            uDisplay[count] = std::move(temp);
        }

        gmsh::view::addModelData(viewTags[4], nbreStep, modelName,
                                 "ElementNodeData", elementTags, uDisplay, t, 3);
    }
}
