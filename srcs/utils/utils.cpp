#include <cassert>
#include "utils.hpp"

bool isPermutation(const std::vector<int>& vec1,
	               const std::vector<int>& vec2,
	               std::vector<unsigned int>& permutation1,
	               std::vector<unsigned int>& permutation2)
{
    assert(vec1.size() != 0 && vec2.size() == vec1.size());
    assert(permutation1.size() == 0 && permutation2.size() == 0);

    permutation1.resize(vec1.size());
    permutation2.resize(vec2.size());

    for(unsigned int i = 0 ; i < vec1.size() ; ++i)
    {
        bool found = false;
        for(unsigned int j = 0 ; j < vec2.size() ; ++j)
        {
            if(vec2[j] == vec1[i])
            {
                //Improve to do not try to check already founded value in vec2.
                permutation1[i] = j;
                permutation2[j] = i;

                found = true;
                break;
            }
        }
        if(!found)
            return false;
    }

    return true;
}


void divideDomain(DomainDiv& domainDiv,
                  int numberProc, const Mesh& mesh)
{
    unsigned int nElements = mesh.nodeData.elementTags.size();

    // Minimum number of nodes per processors
    unsigned int elementDiv = nElements/numberProc;

    // (Possible) remainder of the repartition
    unsigned int remainder = nElements - numberProc * elementDiv;

    // Distribution of the nodes into the processors
    unsigned int offset = 0;

    for(unsigned int i = 0 ; i < numberProc ; i++)
    {
        domainDiv.element[i] = (remainder == 0) ? elementDiv : elementDiv + 1;
        if(remainder != 0) remainder--;

        domainDiv.elementPrec[i] = 0;
        for(unsigned int l = 0 ; l < i ; ++l)
            domainDiv.elementPrec[i] += domainDiv.element[l];

        domainDiv.node[i] = 0;
        domainDiv.nodePrec[i] = 0;

        for(unsigned int k = 0 ; k < domainDiv.element[i] ; ++k)
            domainDiv.node[i]+=mesh.nodeData.elementNumNodes[offset + k];

        offset += domainDiv.element[i]-1;
        for(unsigned j = 0 ; j < i ; ++j)
            domainDiv.nodePrec[i] += domainDiv.node[j];
    }
}
