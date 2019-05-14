#include <gmsh.h>
#include <cmath>
#include "computeNormBis.hpp"
#include "computeNorm.hpp"



void computeNormBis(const Mesh& mesh, const SolverParams& solverParams, double t, Eigen::VectorXd u, double& errorL2, double& errorLinf)
{
    double sum = 0;
    std::vector<double> coeffs = solverParams.initCondition.coefficients;
    std::vector<double> fluxCoeffs = solverParams.fluxCoeffs;

    // loop over the entites
    for(size_t ent = 0 ; ent < mesh.entities.size() ; ++ent)
    {
        // current entity
        Entity entity = mesh.entities[ent];

        //Linf: Vector containing the maximum error among all gauss points, for each element
        Eigen::VectorXd maxGauss(entity.elements.size());
        // loop over the elements
        for(size_t elm = 0 ; elm < entity.elements.size() ; ++elm)
        {
            // current element
            Element element = entity.elements[elm];

            ElementProperty elmProp;

            //get the integration points and basis functions for the same mesh but with a larger number
            //of Gauss points
           gmsh::model::mesh::getBasisFunctions(element.elementTypeHD, "Gauss5",
                                        "Lagrange",
                                        elmProp.intPoints,
                                        elmProp.numComp,
                                        elmProp.basisFunc);

           elmProp.nGP = elmProp.intPoints.size()/4;
           elmProp.nSF = elmProp.basisFunc.size()/elmProp.nGP;

            std::vector<double> PhysIntPointsHD;
            std::vector<double> DeterminantHD;
            std::vector<double> JacobiansHD;

            //Get the integration points in the physical frame and the determinant of the change of variable
            //(x,y,z) => (eta, xi, zeta)
            gmsh::model::mesh::getJacobians(element.elementTypeHD, "Gauss5", JacobiansHD,
                                        DeterminantHD, PhysIntPointsHD, -1);

            //get the coordinates of the current element, for each GP
            std::vector<double> physIntPointsHD(3*elmProp.nGP);
            for (int i = 0; i < 3*elmProp.nGP; ++i)
            {
                    physIntPointsHD[i] = PhysIntPointsHD[3*elmProp.nGP*elm + i];
            }

            //get the determinant of the current element, for each GP
            std::vector<double> determinantHD(elmProp.nGP);
            for (int i = 0; i < elmProp.nGP; ++i)
            {
                    determinantHD[i] = DeterminantHD[elm*elmProp.nGP + i];
            }

            //Linf: Vector containing the error at each Gauss point, for the element elm
            Eigen::VectorXd temp(elmProp.nGP);

            //  loop over the Gauss points
            for(unsigned int k = 0 ; k < elmProp.nGP ; ++k)
            {
                double u_approx = 0;
                double u_exact = 0;
                // loop over the shape functions
                for(unsigned int l = 0 ; l < elmProp.nSF ; ++l)
                {
                    // sum_i{u_i*l_i(x_k)}
                    u_approx += u[element.offsetInU + l]*elmProp.basisFunc[k*elmProp.nSF + l];
                    u_exact += funcGaussian(physIntPointsHD[3*k],t,coeffs,fluxCoeffs)*elmProp.basisFunc[k*elmProp.nSF + l];
                }
                // L2: error = sum_k{w_k*(u(x_k)-u_approx(x_k))^2*det[J](x_k)}
                sum += determinantHD[k]*elmProp.intPoints[4*k+3]*(u_exact-u_approx)*(u_exact-u_approx);

                // L_inf: error  = abs(u(x_k)-u_approx(x_k))
                temp[k] = fabs(u_exact-u_approx);

            }
            // Linf: Maximum over the Gauss points, for the element elm
            maxGauss[elm] = temp.maxCoeff();
        }
        // Linf: Maximum over all elements
        errorLinf = maxGauss.maxCoeff();
        //L2
        errorL2 = sqrt(sum);
    }

}


