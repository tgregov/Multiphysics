#include <cmath>
#include <cassert>
#include "boundaryCondition.hpp"

#include <iostream>

// see .hpp file for description
void reflectShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs)
{

    // check that there is enough values
    assert(field.u.size() == uAtIBC.size());

    // compute a physical reflection
    uAtIBC[0] = field.u[0][indexJ];
    uAtIBC[1] = (-2*edgeNormal[0]*edgeNormal[0]+1)*field.u[1][indexJ]
                    -2*edgeNormal[0]*edgeNormal[1]*field.u[2][indexJ];
    uAtIBC[2] = (-2*edgeNormal[1]*edgeNormal[1]+1)*field.u[2][indexJ]
                    -2*edgeNormal[0]*edgeNormal[1]*field.u[1][indexJ];
}


// see .hpp file for description
void openShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs)
{
    // check that there is enough values
    assert(field.u.size() == uAtIBC.size());

    // compute a physical reflection
    double g = fluxCoeffs[0];
    double H = coeffs[0];
    double alpha = (edgeNormal[0]*field.u[1][indexJ]
                    + edgeNormal[1]*field.u[2][indexJ])/field.u[0][indexJ]
                    - sqrt(g/H)*(field.u[0][indexJ] - H);

    uAtIBC[0] = H;
    uAtIBC[1] = alpha*edgeNormal[0];
    uAtIBC[2] = alpha*edgeNormal[1];
}

// see .hpp file for description
void openAffShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs)
{
    // check that there is enough values
    assert(field.u.size() == uAtIBC.size());

    // compute a physical reflection
    double g = fluxCoeffs[0];
    double H = coeffs[0]*pos[0] + coeffs[1]*pos[1] + coeffs[2];;
    double alpha = (edgeNormal[0]*field.u[1][indexJ]
                        + edgeNormal[1]*field.u[2][indexJ])/field.u[0][indexJ]
                        - sqrt(g/H)*(field.u[0][indexJ] - H);

    uAtIBC[0] = H;
    uAtIBC[1] = alpha*edgeNormal[0];
    uAtIBC[2] = alpha*edgeNormal[1];
}

// see .hpp file for description
void affineShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
              double t, const Field& field, unsigned int indexJ,
              const std::vector<double>& edgeNormal,
              const std::vector<double>& coeffs,
              const std::vector<double>& fluxCoeffs)
{
    uAtIBC[0] = coeffs[0]*pos[0] + coeffs[1]*pos[1] + coeffs[2];
    uAtIBC[1] = 0;
    uAtIBC[2] = 0;
}

// see .hpp file for description
void sinusShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs)
{
    // check that there is enough values
    assert(field.u.size() == uAtIBC.size());
    assert(coeffs.size() == 4);

    // compute a physical reflection
    double g = fluxCoeffs[0];
    double H = coeffs[0]*sin(2*M_PI*coeffs[1]*t + coeffs[2]) + coeffs[3];
    double alpha = (edgeNormal[0]*field.u[1][indexJ]
                    + edgeNormal[1]*field.u[2][indexJ])/field.u[0][indexJ]
                    - sqrt(g/H)*(field.u[0][indexJ] - H);

    uAtIBC[0] = H;
    uAtIBC[1] = alpha*edgeNormal[0];
    uAtIBC[2] = alpha*edgeNormal[1];
}

// see .hpp file for description
void sinusAffShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs)
{
    // check that there is enough values
    assert(field.u.size() == uAtIBC.size());
    assert(coeffs.size() == 6);

    // compute a physical reflection
    double g = fluxCoeffs[0];
    double H = coeffs[0]*sin(2*M_PI*coeffs[1]*t + coeffs[2]) + coeffs[3]
                + coeffs[4]*pos[0] +  coeffs[5]*pos[1];
    double alpha = (edgeNormal[0]*field.u[1][indexJ]
                    + edgeNormal[1]*field.u[2][indexJ])/field.u[0][indexJ]
                    - sqrt(g/H)*(field.u[0][indexJ] - H);

    uAtIBC[0] = H;
    uAtIBC[1] = alpha*edgeNormal[0];
    uAtIBC[2] = alpha*edgeNormal[1];
}

// see .hpp file for description
void gaussian2DShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
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
    uAtIBC[1] = 0;
    uAtIBC[2] = 0;
}


// see .hpp file for description
void gaussian1DShallowX(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs)
{
    // check that there is enough coefficients
    assert(coeffs.size() == 4);
    assert(pos.size() == 3);

    // compute a 1D gaussian along x
    double X = (pos[0]-coeffs[1])*(pos[0]-coeffs[1])/(2*coeffs[2]);

    uAtIBC[0] = coeffs[0]*exp(-(X))+coeffs[3];
    uAtIBC[1] = 0;
    uAtIBC[2] = 0;
}


// see .hpp file for description
void gaussian1DShallowY(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs)
{
    // check that there is enough coefficients
    assert(coeffs.size() == 4);
    assert(pos.size() == 3);

    // compute a 1D gaussian along y
    double Y = (pos[1]-coeffs[1])*(pos[1]-coeffs[1])/(2*coeffs[2]);

    uAtIBC[0] = coeffs[0]*exp(-(Y))+coeffs[3];
    uAtIBC[1] = 0;
    uAtIBC[2] = 0;
}
