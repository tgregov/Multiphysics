#include <gmsh.h>
#include "writer.hpp"

void writeTransport(std::vector<std::vector<double>>& uDisplay,
                  const std::vector<unsigned int>& elementNumNodes,
                  const std::vector<int>& elementTags, const std::string& modelName,
                  unsigned int nbreStep, double t, const Field& field,
                  const std::vector<double>& fluxCoeffs,
                  const std::vector<bool>& whatToWrite, std::vector<int>& viewTags)
{
    if(nbreStep == 0)
    {
        if(whatToWrite[0] == true)
            viewTags[0] = gmsh::view::add("C");
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
}
