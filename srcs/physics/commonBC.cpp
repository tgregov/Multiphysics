#include <cmath>
#include <cassert>
#include "commonBC.hpp"

#include <iostream>


// see .hpp file for description
void constant(std::vector<double>& uAtIBC, const std::vector<double>& pos,
              double t, const Field& field, unsigned int indexJ,
              const std::vector<double>& edgeNormal,
              const std::vector<double>& coeffs,
              const std::vector<double>& fluxCoeffs)
{

    // check that there is enough values
    assert(coeffs.size() == uAtIBC.size());

    // compute constant values
    for(unsigned short unk = 0 ; unk < field.u.size() ; ++unk)
    {
        uAtIBC[unk] = coeffs[unk];
    }
}
