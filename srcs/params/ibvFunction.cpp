#include <cmath>
#include <cassert>
#include "ibvFunction.hpp"

void sinus(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
            const std::vector<double>& u, const std::vector<double>& edgeNormal,
            const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 3);

    uAtIBC[0] = coeffs[0]*sin(2*M_PI*coeffs[1]*t + coeffs[2]);
}

void gaussian(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
                const std::vector<double>& u, const std::vector<double>& edgeNormal,
                const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 3);

    uAtIBC[0] = coeffs[0]*exp(-(t-coeffs[1])*(t-coeffs[1])/(2*coeffs[2]));
}

void constant(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
                const std::vector<double>& u, const std::vector<double>& edgeNormal,
                const std::vector<double>& coeffs)
{
    for(unsigned short unk = 0 ; unk < u.size() ; ++unk)
        uAtIBC[unk] = coeffs[unk];
}

void freeTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
                    const std::vector<double>& u, const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs)
{
    assert(u.size() == uAtIBC.size());

    uAtIBC = u;
}

void reflectShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
                    const std::vector<double>& u, const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs)
{
    assert(u.size() == uAtIBC.size());

    uAtIBC[0] = u[0];
    uAtIBC[1] = (-2*edgeNormal[0]*edgeNormal[0]+1)*u[1]-2*edgeNormal[0]*edgeNormal[1]*u[2];
    uAtIBC[2] = (-2*edgeNormal[1]*edgeNormal[1]+1)*u[2]-2*edgeNormal[0]*edgeNormal[1]*u[1];
}

void gaussian2DShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
                        const std::vector<double>& u, const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 6);
    assert(pos.size() == 3);

    double X = (pos[0]-coeffs[1])*(pos[0]-coeffs[1])/(2*coeffs[2]);
    double Y = (pos[1]-coeffs[3])*(pos[1]-coeffs[3])/(2*coeffs[4]);

    uAtIBC[0] = coeffs[0]*exp(-(X+Y))+coeffs[5];
    uAtIBC[1] = 0;
    uAtIBC[2] = 0;
}

void gaussian1DShallowX(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
                        const std::vector<double>& u, const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 4);
    assert(pos.size() == 3);

    double X = (pos[0]-coeffs[1])*(pos[0]-coeffs[1])/(2*coeffs[2]);

    uAtIBC[0] = coeffs[0]*exp(-(X))+coeffs[3];
    uAtIBC[1] = 0;
    uAtIBC[2] = 0;
}

void gaussian1DShallowY(std::vector<double>& uAtIBC, const std::vector<double>& pos, double t,
                        const std::vector<double>& u, const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 4);
    assert(pos.size() == 3);

    double Y = (pos[1]-coeffs[1])*(pos[1]-coeffs[1])/(2*coeffs[2]);

    uAtIBC[0] = coeffs[0]*exp(-(Y))+coeffs[3];
    uAtIBC[1] = 0;
    uAtIBC[2] = 0;
}
