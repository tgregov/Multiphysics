#include <gmsh.h>
#include "write.hpp"

void writeEnd(const std::vector<int>& viewTags, const std::vector<bool>& whatToWrite,
                const std::string& resultsName)
{
    for(unsigned int i = 0 ; i < whatToWrite.size() ; ++i)
    {
        if(whatToWrite[i] == true)
            gmsh::view::write(viewTags[i], resultsName, true);
    }
}
