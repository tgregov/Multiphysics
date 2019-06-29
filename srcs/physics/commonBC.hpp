#ifndef commonBC_hpp_included
#define commonBC_hpp_included

#include <vector>
#include "../solver/field.hpp"

/**
 * \brief Compute a constant -- for shallow waters & the pure transport case.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficient for the constant: coeffs[unk] = constant.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void constant(std::vector<double>& uAtIBC, const std::vector<double>& pos,
              double t, const Field& field, unsigned int indexJ,
              const std::vector<double>& edgeNormal,
              const std::vector<double>& coeffs,
              const std::vector<double>& fluxCoeffs);

#endif // commonBC_hpp_included
