#include <cmath>
#include <cassert>
#include "ibvFunction.hpp"

double sinus(const std::vector<double>& pos, double u,
             double t, const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 3);

    return coeffs[0]*sin(2*M_PI*coeffs[1]*t + coeffs[2]);
}

double gaussian(const std::vector<double>& pos, double u,
                double t, const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 3);

    return coeffs[0]*exp(-(t-coeffs[1])*(t-coeffs[1])/(2*coeffs[2]));
}

double constant(const std::vector<double>& pos, double u,
                double t, const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 1);

    return coeffs[0];
}

double constantNeumann(const std::vector<double>& pos, double u,
                double t, const std::vector<double>& coeffs)
{
    return u;
}

double gaussian2D(const std::vector<double>& pos, double u,
                  double t, const std::vector<double>& coeffs)
{
    assert(coeffs.size() == 5);
    assert(pos.size() == 3);

    double X = (pos[0]-coeffs[1])*(pos[0]-coeffs[1])/(2*coeffs[2]);
    double Y = (pos[1]-coeffs[3])*(pos[1]-coeffs[3])/(2*coeffs[4]);

    return coeffs[0]*exp(-(X+Y));
}
