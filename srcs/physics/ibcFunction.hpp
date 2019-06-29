#ifndef ibcFunction_hpp_included
#define ibcFunction_hpp_included

#include <functional>
#include <vector>
#include "../solver/field.hpp"

/**
 * \struct bc
 * \brief Mathematical function for a boundary condition.
 */
struct ibc
{
    std::vector<double> coefficients; /**< Coefficient for the mathematical function */
    std::function<void(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs)> ibcFunc;
                        /**< Pointer to the mathematical function */
};

#endif // ibcFunction_hpp_included

