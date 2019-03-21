#include <cmath>
#include <cassert>
#include "bcFunction.hpp"

double sinus(double x, double y, double z,
             double t, std::vector<double> coeffs)
{
    assert(coeffs.size() == 3);

    return coeffs[0]*sin(coeffs[1]*t + coeffs[2]);
}

double gaussian(double x, double y, double z,
                double t, std::vector<double> coeffs)
{
    assert(coeffs.size() == 3);

    return coeffs[0]*exp(-(t-coeffs[1])*(t-coeffs[1])/coeffs[2]);
}

double constant(double x, double y, double z,
                double t, std::vector<double> coeffs)
{
    assert(coeffs.size() == 1);

    return coeffs[0];
}
