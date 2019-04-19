#include <cassert>
#include <vector>
#include <iostream>

bool isPermutation(const std::vector<int>& vec1,
	               const std::vector<int>& vec2,
	               std::vector<unsigned int>& permutation1,
	               std::vector<unsigned int>& permutation2);

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
