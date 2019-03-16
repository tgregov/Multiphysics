#include <iostream>
#include "buildS.hpp"

void buildS(const Mesh2D& mesh, Eigen::SparseMatrix<double>& Sx,
			Eigen::SparseMatrix<double>& Sy)
{

    // indexx, indexy are vectors of triplets that contains the coordinates in the
    // [Sx], [Sy] matrices of each ot their components; the offsetMatrix is the
    // upper-left coordinate at which the current element matrix should be added
	std::vector<Eigen::Triplet<double>> indexx, indexy;
    unsigned int offsetMatrix = 0;

    // loop over the entites
    for(unsigned int ent = 0 ; ent < mesh.entities.size() ; ++ent)
    {
        Entity2D entity = mesh.entities[ent];

        // loop over the elements
        for(unsigned int elm = 0 ; elm < entity.elements.size() ; ++elm)
        {
            Element2D element = entity.elements[elm];

            // get the properties of the current element type
            ElementProperty elmProp = mesh.elementProperties2D.at(element.elementType2D);
            std::vector<std::vector<double>> pondFunc = elmProp.pondFunc;

            // matrices [Sx], [Sy] for the current element, represented as vectors
            std::vector<double> sxElm(elmProp.nSF*elmProp.nSF, 0.0);
			std::vector<double> syElm(elmProp.nSF*elmProp.nSF, 0.0);

			for(unsigned int k = 0 ; k < elmProp.nGP ; ++k)
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
				double dxdxi = element.jacobian2D[9*k];
				double dxdeta = element.jacobian2D[9*k + 3];
				double dydxi = element.jacobian2D[9*k + 1];
				double dydeta = element.jacobian2D[9*k + 4];

				// N.B.: the dets simplify in the inverse of a 2x2 matrix !
				double dxidx = dydeta;
				double detadx = - dxdeta;
				double dxidy = - dydxi;
				double detady = dxdxi;

				for(unsigned int i = 0; i < elmProp.nSF; i++)
				{
					for(unsigned int j = 0; j < elmProp.nSF; j++)
					{

						double dljdxi = elmProp.basisFuncGrad
							[k*elmProp.nSF*3 + j*3];
						double dljdeta = elmProp.basisFuncGrad
							[k*elmProp.nSF*3 + j*3 + 1];

						sxElm[i*elmProp.nSF + j]
							+= pondFunc[k][i]*(dljdxi*dxidx + dljdeta*detadx);
						syElm[i*elmProp.nSF + j]
							+= pondFunc[k][i]*(dljdxi*dxidy + dljdeta*detady);


						// if we have calculated the sum for all the GP, we can
						// save the computed components
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

			offsetMatrix += elmProp.nSF;
		}
   	}

    Sx.setFromTriplets(indexx.begin(), indexx.end());
    Sy.setFromTriplets(indexy.begin(), indexy.end());
}
