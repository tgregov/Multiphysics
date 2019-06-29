#include <iostream>
#include "buildS.hpp"


// see .hpp for description
void buildS(const Mesh& mesh, Eigen::SparseMatrix<double>& Sx,
			Eigen::SparseMatrix<double>& Sy)
{

	// * indexx/indexy: vectors of triplets that contains the coordinates in the
    // 		[Sx]/[Sy] matrices of each ot their components
    // * offsetMatrix: upper-left coordinate at which the current element matrix
    //		should be added
	std::vector<Eigen::Triplet<double>> indexx, indexy;
    unsigned int offsetMatrix = 0;

    // loop over the elements
    for(size_t elm = 0 ; elm < mesh.elements.size() ; ++elm)
    {
        // get the 2D properties of the current element type
        // * pondFunc[k][i]: w_k*l_i evaluated at each GP
        ElementProperty elmProp
            = mesh.elementProperties.at(mesh.elements[elm].elementTypeHD);
        std::vector<std::vector<double>> pondFunc = elmProp.pondFunc;

        // matrices [Sx], [Sy] for the current element, the matrices are stored
        // as vectors, such that S_{i,j} = [S](i*elmProp.nSF + j)
        std::vector<double> sxElm(elmProp.nSF*elmProp.nSF, 0.0);
        std::vector<double> syElm(elmProp.nSF*elmProp.nSF, 0.0);

        // sum over the GP
        for(unsigned int k = 0 ; k < elmProp.nGP ; ++k)
        {

            // only the 2D case here
            double dxdxi 	= mesh.elements[elm].jacobianHD[9*k];
            double dxdeta 	= mesh.elements[elm].jacobianHD[9*k + 3];
            double dydxi 	= mesh.elements[elm].jacobianHD[9*k + 1];
            double dydeta 	= mesh.elements[elm].jacobianHD[9*k + 4];

            // dzdzeta component of the jacobian: +/-1
            double sign = mesh.elements[elm].jacobianHD[9*k + 8];

            // the dets simplify in the inverse of a 2x2 matrix !
            double dxidx 	= dydeta/sign;
            double detadx 	= - dydxi/sign;
            double dxidy 	= - dxdeta/sign;
            double detady 	= dxdxi/sign;

            // loop over the components of the matrix
            for(unsigned int i = 0 ; i < elmProp.nSF ; ++i)
            {
                for(unsigned int j = 0 ; j < elmProp.nSF ; ++j)
                {

                    // gradient of the shape functions
                    double dljdxi = elmProp.basisFuncGrad
                        [k*elmProp.nSF*3 + j*3];
                    double dljdeta = elmProp.basisFuncGrad
                        [k*elmProp.nSF*3 + j*3 + 1];

                    // components of the elements
                    sxElm[i*elmProp.nSF + j]
                        += pondFunc[k][i]*(dljdxi*dxidx + dljdeta*detadx);
                    syElm[i*elmProp.nSF + j]
                        += pondFunc[k][i]*(dljdxi*dxidy + dljdeta*detady);

                    // if we have calculated the sum for all the GP,
                    // we can save the computed components
                    if(k == elmProp.nGP - 1)
                    {
                        indexx.push_back(Eigen::Triplet<double>
                            (i + offsetMatrix, j + offsetMatrix,
                                sxElm[i*elmProp.nSF + j]));

                        indexy.push_back(Eigen::Triplet<double>
                            (i + offsetMatrix, j + offsetMatrix,
                                syElm[i*elmProp.nSF + j]));
                    }
                }
            }
        }

        // increase the offset of the local matrix
        offsetMatrix += elmProp.nSF;
    }

    // add the triplets in the sparse matrix
    Sx.setFromTriplets(indexx.begin(), indexx.end());
    Sy.setFromTriplets(indexy.begin(), indexy.end());
}
