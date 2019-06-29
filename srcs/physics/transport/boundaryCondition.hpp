#ifndef transport_boundaryCondition_hpp_included
#define transport_boundaryCondition_hpp_included

#include <vector>
#include "../../solver/field.hpp"

/**
 * \brief Compute a wave of the shape A*sin(2*pi*nu*t + phi) -- for the pure
 * transport case.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the sinus: coeffs[0] = A, coeffs[1] = nu,
 * coeffs[2] = phi.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void sinusTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
           double t, const Field& field, unsigned int indexJ,
           const std::vector<double>& edgeNormal,
           const std::vector<double>& coeffs,
           const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a gaussian: A*exp(-(t-t_peak)^2/(2*var)) -- for the pure transport
 * case.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = t_peak,
 * coeffs[2] = var.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussianTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
              double t, const Field& field, unsigned int indexJ,
              const std::vector<double>& edgeNormal,
              const std::vector<double>& coeffs,
              const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a 2D gaussian:  A*exp(-(x-x0)^2/(2*var_y)-(y-y0)^2/(2*var_y)) --
 * for the pure transport case.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = y_0, coeffs[4] = var_y.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussian2DTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                         double t, const Field& field, unsigned int indexJ,
                         const std::vector<double>& edgeNormal,
                         const std::vector<double>& coeffs,
                         const std::vector<double>& fluxCoeffs);


/**
 * \brief Compute a physical opening -- for the pure transport case.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficient (not used here).
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void freeTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                   double t, const Field& field, unsigned int indexJ,
                   const std::vector<double>& edgeNormal,
                   const std::vector<double>& coeffs,
                   const std::vector<double>& fluxCoeffs);

#endif // transport_boundaryCondition_hpp_included
