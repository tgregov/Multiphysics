#ifndef shallow_boundaryCondition_hpp_included
#define shallow_boundaryCondition_hpp_included

#include <vector>
#include "../../solver/field.hpp"



/**
 * \brief Compute a wave of the shape A*sin(2*pi*nu*t + phi) + B -- for the
 * shallow water.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the sinus: coeffs[0] = A, coeffs[1] = nu,
 * coeffs[2] = phi, coeffs[3] = B.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void sinusShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a physical opening -- for shallow waters.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient, coeffs[0] = H away from the BC.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void openAffShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a linear initial condtion -- for shallow waters.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient, coeffs[0] = dh/dx, coeffs[0] = dh/dy, coeffs[0] = offset.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void affineShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
              double t, const Field& field, unsigned int indexJ,
              const std::vector<double>& edgeNormal,
              const std::vector<double>& coeffs,
              const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a sinus BC when there is a seabed gradient -- for shallow waters.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient, coeffs[0] = A , coeffs[1] = f , coeffs[2] = phi,
 coeffs[3] = offset, coeffs[] = dh/dx, coeffs[0] = dh/dy
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void sinusAffShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a physical reflection -- for shallow waters.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficient (not used here).
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void reflectShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a physical opening -- for shallow waters.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient (not used here).
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void openShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs,
                    const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a 2D gaussian: A*exp(-(x-x0)^2/(2*var_y)-(y-y0)^2/(2*var_y)) -- for
 * shallow waters
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = y_0, coeffs[4] = var_y.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussian2DShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs);

/**
 * \brief Compute a 1D gaussian along x: A*exp(-(x-x0)^2/(2*var_x) + B -- for shallow
 * waters.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = B.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussian1DShallowX(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs);


/**
 * \brief Compute a 1D gaussian along y: A*exp(-(y-y0)^2/(2*var_y) + B -- for shallow
 * waters.
 * \param uAtIBC Value of the BC unknowns at (x, y, z, t).
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = y_0,
 * coeffs[2] = var_y, coeffs[3] = B.
 * \param fluxCoeffs Coefficients of the physical fluxes.
 */
void gaussian1DShallowY(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs,
                        const std::vector<double>& fluxCoeffs);


#endif // shallow_boundaryCondition_hpp_included
