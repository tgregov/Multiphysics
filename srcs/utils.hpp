#ifndef permutation_hpp_included
#define permutation_hpp_included

#include <vector>

bool isPermutation(const std::vector<int>& vec1,
	               const std::vector<int>& vec2,
	               std::vector<unsigned int>& permutation1,
	               std::vector<unsigned int>& permutation2);

#endif /* permutation_hpp_included */