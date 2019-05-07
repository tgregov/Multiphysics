#ifndef utils_hpp_included
#define utils_hpp_included

#include <vector>
#include "../mesh/Mesh.hpp"

struct DomainDiv
{
    std::vector<int> element;
    std::vector<int> elementPrec;
    std::vector<int> node;
    std::vector<int> nodePrec;

    DomainDiv(unsigned int numberProc)
    {
        element.resize(numberProc);
        elementPrec.resize(numberProc);
        node.resize(numberProc);
        nodePrec.resize(numberProc);
    }
};

bool isPermutation(const std::vector<int>& vec1,
	               const std::vector<int>& vec2,
	               std::vector<unsigned int>& permutation1,
	               std::vector<unsigned int>& permutation2);

void divideDomain(DomainDiv& domainDiv,
                  int numberProc, const Mesh& mesh);

#endif /* utils_hpp_included */
