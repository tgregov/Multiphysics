#include <gmsh.h>
#include "computeL2Norm.hpp"

double func(double x, double t, const std::vector<double>& coeffs, const std::vector<double>& fluxCoeffs){

    double v = sqrt(fluxCoeffs[0]*fluxCoeffs[1]);

return coeffs[3] + 0.5*(coeffs[0]*exp(-(x-coeffs[1]+v*t)*(x-coeffs[1]+v*t)/(2*coeffs[2]))+ 
    coeffs[0]*exp(-(x-coeffs[1]-v*t)*(x-coeffs[1]-v*t)/(2*coeffs[2])));}


double computeL2Norm(const Mesh& mesh, const SolverParams& solverParams, double t, Eigen::VectorXd u)
{
    double sum = 0;

    std::vector<double> coeffs = solverParams.initCondition.coefficients;
    std::vector<double> fluxCoeffs = solverParams.fluxCoeffs;

    // loop over the entites
    for(size_t ent = 0 ; ent < mesh.entities.size() ; ++ent)
    {
        // current entity
        Entity entity = mesh.entities[ent];

        // loop over the elements
        for(size_t elm = 0 ; elm < entity.elements.size() ; ++elm)
        {
            // current element
            Element element = entity.elements[elm];

            ElementProperty elmProp;

            //get the integration points and basis functions for the samemesh but with a given number
            //of Gauss points
           gmsh::model::mesh::getBasisFunctions(element.elementTypeHD, "Gauss30",
                                        "Lagrange",
                                        elmProp.intPoints,
                                        elmProp.numComp,
                                        elmProp.basisFunc);
           //std::cout << elmProp.intPoints.size() << std::endl;

           elmProp.nGP = elmProp.intPoints.size()/4;
           elmProp.nSF = elmProp.basisFunc.size()/elmProp.nGP;

            std::vector<double> PhysIntPointsHD;
            std::vector<double> DeterminantHD;
            std::vector<double> JacobiansHD;

            //Get the integration points in the physical frame and the determinant of the change of variable
            //(x,y,z) => (eta, xi, zeta)
            gmsh::model::mesh::getJacobians(element.elementTypeHD, "Gauss30", JacobiansHD,
                                        DeterminantHD, PhysIntPointsHD, -1);

            //get the coordinates of the current element
            std::vector<double> physIntPointsHD(3*elmProp.nGP);
            for (int i = 0; i < 3*elmProp.nGP; ++i)
            {
                    physIntPointsHD[i] = PhysIntPointsHD[3*elmProp.nGP*elm + i];
            }

            //get the determinant of the current element
            std::vector<double> determinantHD(elmProp.nGP);
            for (int i = 0; i < elmProp.nGP; ++i)
            {
                    determinantHD[i] = DeterminantHD[elm*elmProp.nGP + i];
            }


            //  loop over the Gauss points
            for(unsigned int k = 0 ; k < elmProp.nGP ; ++k)
            {
                double u_approx = 0;
                // loop over the shape functions
                for(unsigned int l = 0 ; l < elmProp.nSF ; ++l)
                {
                    // sum_i{u_i*l_i(x_k)}
                    u_approx += u[element.offsetInU + l]*elmProp.basisFunc[k*elmProp.nSF + l];
                }
                // error = sum_k{w_k*(u(x_k)-u_approx(x_k))^2*det[J](x_k)}
                sum += determinantHD[k]*elmProp.intPoints[4*k+3]
                    *(func(physIntPointsHD[3*k],t,coeffs,fluxCoeffs)-u_approx)
                    *(func(physIntPointsHD[3*k],t,coeffs,fluxCoeffs)-u_approx);
            }
        }
    }

    return sum;
}