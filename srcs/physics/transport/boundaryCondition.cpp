#include <cmath>
#include <cassert>
#include "boundaryCondition.hpp"

#include <iostream>

// see .hpp file for description
void sinusTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
           double t, const Field& field, unsigned int indexJ,
           const std::vector<double>& edgeNormal,
           const std::vector<double>& coeffs,
           const std::vector<double>& fluxCoeffs)
{
    // check that there is enough coefficients
    assert(coeffs.size() == 4);

    // compute a sine wave
    uAtIBC[0] = coeffs[0]*sin(2*M_PI*coeffs[1]*t + coeffs[2])+coeffs[3];
}


// see .hpp file for description
void gaussianTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
           double t, const Field& field, unsigned int indexJ,
           const std::vector<double>& edgeNormal,
           const std::vector<double>& coeffs,
           const std::vector<double>& fluxCoeffs)
{
    // check that there is enough coefficients
    assert(coeffs.size() == 4);

    // compute a sine wave
    uAtIBC[0] = coeffs[0]*sin(2*M_PI*coeffs[1]*t + coeffs[2])+coeffs[3];
}


// see .hpp file for description
void freeTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                   double t, const Field& field, unsigned int indexJ,
                   const std::vector<double>& edgeNormal,
                   const std::vector<double>& coeffs,
                   const std::vector<double>& fluxCoeffs)
{

    // check that there is enough values
    assert(field.u.size() == uAtIBC.size());

    // compute same values
    for(unsigned short unk = 0 ; unk < field.u.size() ; ++unk)
        uAtIBC[unk] = field.u[unk][indexJ];
}


// see .hpp file for description
void gaussian2DTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                         const std::vector<double>& fluxCoeffs)
{
    // check that there is enough coefficients
    assert(coeffs.size() == 6);
    assert(pos.size() == 3);

    // compute a 2D gaussian
    double X = (pos[0]-coeffs[1])*(pos[0]-coeffs[1])/(2*coeffs[2]);
    double Y = (pos[1]-coeffs[3])*(pos[1]-coeffs[3])/(2*coeffs[4]);

    uAtIBC[0] = coeffs[0]*exp(-(X+Y))+coeffs[5];
}
