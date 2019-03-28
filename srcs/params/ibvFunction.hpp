#ifndef bcFunction_hpp_included
#define bcFunction_hpp_included

#include <functional>
#include <map>
#include <vector>

/**
 * \struct bc
 * \brief Mathematical function for a boundary condition.
 */
struct ibc
{
    std::vector<double> coefficients;   /**< Coefficient for the mathematical function */
    std::function<double(double x, double y, double z,
                         double u, double t,
                         std::vector<double> coeffs)> ibcFunc; /**< Pointer to the mathematical function */
};


/**
 * \brief Compute a simple A*sin(2*pi*nu*t + phi)
 * \param x x coordinate.
 * \param y y coordinate.
 * \param z z coordinate.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient for the sinus. coeffs[0] = A, coeffs[1] = nu
 * coeffs[2] = phi.
 * \return Value of the function at (x, y , z, t).
 */
double sinus(double x, double y, double z, double u,
             double t, const std::vector<double>& coeffs);

/**
 * \brief Compute a simple A*exp(-(t-t_peak)^2/(2*var))
 * \param x x coordinate.
 * \param y y coordinate.
 * \param z z coordinate.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient for the gaussian. coeffs[0] = A, coeffs[1] = t_peak
 * coeffs[2] = var.
 * \return Value of the function at (x, y , z, t).
 */
double gaussian(double x, double y, double z, double u,
                double t, const std::vector<double>& coeffs);

/**
 * \brief Compute a constant
 * \param x x coordinate.
 * \param y y coordinate.
 * \param z z coordinate.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient for the constant. coeffs[0] = constant.
 * \return Value of the function at (x, y , z, t).
 */
double constant(double x, double y, double z, double u,
                double t, const std::vector<double>& coeffs);

/**
 * \brief Compute a Von Neumann constant value
 * \param x x coordinate.
 * \param y y coordinate.
 * \param z z coordinate.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient (not used here).
 * \return u.
 */
double constantNeumann(double x, double y, double z, double u,
                       double t, const std::vector<double>& coeffs);

/**
 * \brief Compute a simple A*exp(-(x-x0)^2/(2*var_y) -(y-y0)^2/(2*var_y))
 * \param x x coordinate.
 * \param y y coordinate.
 * \param z z coordinate.
 * \param u The current solution.
 * \param t Current time.
 * \param coeffs Coefficient for the gaussian. coeffs[0] = A, coeffs[1] = x_0
 * coeffs[2] = var_x, coeffs[3]=y_0, coeffs[4]=var_y.
 * \return Value of the function at (x, y , z, t).
 */
double gaussian2D(double x, double y, double z, double u,
                  double t, const std::vector<double>& coeffs);

#endif // bcFunction_hpp_included
