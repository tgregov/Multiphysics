#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
#include <Eigen/Sparse>
#include "buildS.hpp"

/* Function that builds the [Sx] and [Sy] matrix of the DG method
 * Inputs:
 *  - meshParams: structure that contains all the geometric & mesh information
 *  - Sx: sparse matrix [Sx], defined before this function call
 *  - Sy: sparse matrix [Sx], defined before this function call
 *
 * Description:
 * The objective is to compute [Sx] and [Sy], defined respectively by
 * sum_k{w_k*l_i(x_k)*[dl_j/dxi(x_k)*dxi/dx(x_k) 
 *						+ dl_j/deta(x_k)*deta/dx(x_k)]*det[J](x_k)}, and
 $ sum_k{w_k*l_i(x_k)*[dl_j/dxi(x_k)*dxi/dy(x_k) 
 *						+ dl_j/deta(x_k)*deta/dy(x_k)]*det[J](x_k)}, where the
 * sum is done over the Gauss points (GP). This is done in 3 steps:
 *  I. Precomputation of what is element-independent:
 *		More precisely, the following components are considered: w_k*l_i(x_k)
 *		are stored in a 2D vector WL such that WL[k][i] = w_k*l_i(x_k). The 
 *		derivatives of the shape functions, that is, dl_j/dxi and dl_j/deta are
 *		already contained in meshParams.basisFuncGrad.
 *  II. Computation of what is element-dependent:
 *      We define by sx/sy to be vector of vectors, where sx/sy are vectors (of 
 *		vectors) of length nE, and sx[0]/sy[0] vectors of length nSF*nSF, each 
 *		component being concatenated (first loop over i, then over j). 
 *		sxElm/syElm is defined as temporary vector such that when element elm is
 *		considered, sxElm/syElm equals sx[elm]/sy[elm]. 
 *		These elements are computed in the following way: a first loop is done
 *		over the elements, so that we consider the computation of sxElm/syElm.
 *		Then, a loop is done over the GP, to compute an element relative to a
 *		a specific k in the original formula. For that GP, the inverse of the
 *		jacobian matrix is determined. Then, still for that GP, all the (i,j)
 *		components are computed. When looping over the GP, the intermediate 
 *		result k is added to all the previous ones, to get in the end the sum 
 *		over all the GP, i.e. over k. Finally, the vectors sxElm/syElm is added
 *		to the vectors of vectors sx/sy.
 *  III. Storage of the result using Eigen:
 *      The matrices [Sx]/[Sy] are assembled using the data contained in sx/sy. 
 *		The order is determined using the fact that sx/sy are vectors such that
 *		sx[elm]/sy[elm] contains the (i,j) concanated components relative to 
 *		the element elm.
 */     

void buildS(const MeshParams& meshParams, Eigen::SparseMatrix<double>& Sx, 
			Eigen::SparseMatrix<double>& Sy)
{

    // I. Precomputation of what is element-independent
	std::vector<std::vector<double>> WL;

	for(unsigned int k = 0; k < meshParams.nGP; k++)
	{
		std::vector<double> wl;

		for(unsigned int i = 0; i < meshParams.nSF; i++)
		{
			wl.push_back(meshParams.intPoints[4*k + 3]
							*meshParams.basisFunc[meshParams.nSF*k + i]);
			// meshParams.intPoints[4*k + 3] = w_k
			// meshParams.basisFunc[meshParams.nSF*k + i] = l_i(x_k)
		}

		WL.push_back(wl);
	}

	// II. Computation of what is element-dependent
	std::vector<std::vector<double>> sx, sy;

	for(unsigned int elm = 0; elm < meshParams.nE; elm++)
	{

		std::vector<double> sxElm(meshParams.nSF*meshParams.nSF, 0.0);
		std::vector<double> syElm(meshParams.nSF*meshParams.nSF, 0.0);

		for(unsigned int k = 0; k < meshParams.nGP; k++)
		{	
			// only the 2D case here; ideally it should be more general
			// [TO CHECK] in the SDK, it is said that jacobian is stored as
			// [Jxx, Jxy, Jxz, ...]. However, it seems that it is reversed:
			// by manually checking, I get that jacobian[0] = dx1/dX1,
			// jacobian[1] = dx2/dX1, etc
			// N.B.1: maybe it depends on how we define the jacobian matrix [J]
			// N.B.2: how to check manually ? We have our simple 4 triangles
			// in a square geometry. We can assume a linear change of variables
			// => 	x = (dx/dxi)*xi + (dx/deta)*eta + (dx/dzeta)*zeta + kappa_x
			// 		y = (dy/dxi)*xi + (dy/deta)*eta + (dy/dzeta)*zeta + kappa_y
			//		z = (dz/dxi)*xi + (dz/deta)*eta + (dz/dzeta)*zeta + kappa_z
			// then, just take (xi, eta) in {(0, 0), (1, 0), (0, 1)} and verify
			// that it works (the kappas correspond to simple translation)
			// N.B.3: the basisFuncGrad seems alright by checking manually
			double dxdxi = meshParams.jacobian[9*meshParams.nGP*elm + 9*k];
			double dxdeta = meshParams.jacobian[9*meshParams.nGP*elm + 9*k + 3];
			double dydxi = meshParams.jacobian[9*meshParams.nGP*elm + 9*k + 1];
			double dydeta = meshParams.jacobian[9*meshParams.nGP*elm + 9*k + 4];

			// N.B.: the dets simplify in the inverse of a 2x2 matrix !
			double dxidx = dydeta;
			double detadx = - dxdeta;
			double dxidy = - dydxi;
			double detady = dxdxi;

			for(unsigned int i = 0; i < meshParams.nSF; i++)
			{
				for(unsigned int j = 0; j < meshParams.nSF; j++)
				{

					double dljdxi = meshParams.basisFuncGrad
						[k*meshParams.nSF*3 + j*3];
					double dljdeta = meshParams.basisFuncGrad
						[k*meshParams.nSF*3 + j*3 + 1];

					sxElm[i* meshParams.nSF + j] 
						+= WL[k][i]*(dljdxi*dxidx + dljdeta*detadx);
					syElm[i* meshParams.nSF + j] 
						+= WL[k][i]*(dljdxi*dxidy + dljdeta*detady);
				}
			}
		}

		sx.push_back(sxElm);
		sy.push_back(syElm);
	}

	// III. Storage of the result using Eigen
    std::vector<Eigen::Triplet<double>> indexx, indexy;

    for(unsigned int elm = 0; elm < meshParams.nE; elm++)
    {
    	for(unsigned int i = 0; i < meshParams.nSF; i++)
    	{
    		for(unsigned int j = 0; j < meshParams.nSF; j++)
    		{
            	indexx.push_back(Eigen::Triplet<double>
                	(i + elm*meshParams.nSF, j+ elm*meshParams.nSF, 
                		sx[elm][i*meshParams.nSF + j]));

    			indexy.push_back(Eigen::Triplet<double>
                	(i + elm*meshParams.nSF, j+ elm*meshParams.nSF, 
                		sy[elm][i*meshParams.nSF + j]));
    		}
    	}
    }

    Sx.setFromTriplets(indexx.begin(), indexx.end());
    Sy.setFromTriplets(indexy.begin(), indexy.end());
}