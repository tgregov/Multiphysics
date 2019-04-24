#ifndef bcFunction_hpp_included
#define bcFunction_hpp_included

#include <functional>
#include <map>
#include <vector>
#include "../solver/field.hpp"


/**
 * \struct bc
 * \brief Mathematical function for a boundary condition.
 */
struct ibc
{
    std::vector<double> coefficients;   /**< Coefficient for the mathematical function */
    std::function<void(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs)> ibcFunc; /**< Pointer to the mathematical function */
};


/**
 * \brief Compute a wave of the shape A*sin(2*pi*nu*t + phi) -- for the pure
 * transport case.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the sinus: coeffs[0] = A, coeffs[1] = nu,
 * coeffs[2] = phi.
 * \return Value of the function at (x, y, z, t).
 */
void sinus(std::vector<double>& uAtIBC, const std::vector<double>& pos,
           double t, const Field& field, unsigned int indexJ,
           const std::vector<double>& edgeNormal,
           const std::vector<double>& coeffs);


/**
 * \brief Compute a gaussian: A*exp(-(t-t_peak)^2/(2*var)) -- for the pure transport
 * case.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = t_peak,
 * coeffs[2] = var.
 * \return Value of the function at (x, y, z, t).
 */
void gaussian(std::vector<double>& uAtIBC, const std::vector<double>& pos,
              double t, const Field& field, unsigned int indexJ,
              const std::vector<double>& edgeNormal,
              const std::vector<double>& coeffs);


/**
 * \brief Compute a 2D gaussian:  A*exp(-(x-x0)^2/(2*var_y)-(y-y0)^2/(2*var_y)) --
 * for the pure transport case.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = y_0, coeffs[4] = var_y.
 * \return Value of the function at (x, y, z, t).
 */
void gaussian2DTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs);


/**
 * \brief Compute a Neumann constant value -- for the pure transport case.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficient (not used here).
 * \return u.
 */
void freeTransport(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                   double t, const Field& field, unsigned int indexJ,
                   const std::vector<double>& edgeNormal,
                   const std::vector<double>& coeffs);


/**
 * \brief Compute a constant -- for shallow waters & the pure transport case.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficient for the constant: coeffs[0] = constant.
 * \return Value of the function at (x, y, z, t).
 */
void constant(std::vector<double>& uAtIBC, const std::vector<double>& pos,
              double t, const Field& field, unsigned int indexJ,
              const std::vector<double>& edgeNormal,
              const std::vector<double>& coeffs);


/**
 * \brief Compute a physical reflection -- for shallow waters.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficient (not used here).
 * \return Physically reflected values.
 */
void reflectShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                    double t, const Field& field, unsigned int indexJ,
                    const std::vector<double>& edgeNormal,
                    const std::vector<double>& coeffs);


/**
 * \brief Compute a 2D gaussian: A*exp(-(x-x0)^2/(2*var_y)-(y-y0)^2/(2*var_y)) -- for
 * shallow waters
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = y_0, coeffs[4] = var_y.
 * \return Value of the function at (x, y, z, t).
 */
void gaussian2DShallow(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs);


/**
 * \brief Compute a 1D gaussian along x: A*exp(-(x-x0)^2/(2*var_x) + B -- for shallow
 * waters.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = x_0,
 * coeffs[2] = var_x, coeffs[3] = B.
 * \return Value of the function at (x, y, z, t).
 */
void gaussian1DShallowX(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs);


/**
 * \brief Compute a 1D gaussian along y: A*exp(-(y-y0)^2/(2*var_y) + B -- for shallow
 * waters.
 * \param pos Node position.
 * \param t Current time.
 * \param field Structure containing the current solution.
 * \param indexJ Index of the node corresponding to the boundary in the field structure.
 * \param coeffs Coefficients for the gaussian: coeffs[0] = A, coeffs[1] = y_0,
 * coeffs[2] = var_y, coeffs[3] = B.
 * \return Value of the function at (x, y, z, t).
 */
void gaussian1DShallowY(std::vector<double>& uAtIBC, const std::vector<double>& pos,
                        double t, const Field& field, unsigned int indexJ,
                        const std::vector<double>& edgeNormal,
                        const std::vector<double>& coeffs);

#endif // bcFunction_hpp_included
