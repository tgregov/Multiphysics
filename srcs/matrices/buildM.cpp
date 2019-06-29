#include <iostream>
#include "buildM.hpp"


// see .hpp file for description
void buildM(const Mesh& mesh, Eigen::SparseMatrix<double>& invM)
{

    // * index: vector of triplets that contains the coordinates in the [M] matrix
    //          for each of its component
    // * offsetMatrix: upper-left coordinate at which the current element matrix
    //          should be added
    std::vector<Eigen::Triplet<double>> index;
    unsigned int offsetMatrix = 0;

    // loop over the elements
    for(size_t elm = 0 ; elm < mesh.elements.size() ; ++elm)
    {
        // get the 2D properties of the current element type
        // * prodFunc[k][i,j]: w_k*l_i*l_j evaluated at each GP
        // * IJ[l]: components (i, j) for the index l that runs through the
        //          upper-half part of [M]
        ElementProperty elmProp
            = mesh.elementProperties.at(mesh.elements[elm].elementTypeHD);
        std::vector<std::vector<double>> prodFunc = elmProp.prodFunc;
        std::vector<std::pair<unsigned int, unsigned int>> IJ = elmProp.IJ;

        // local [M] matrix for the current element
        Eigen::MatrixXd MLocal(elmProp.nSF, elmProp.nSF);
        MLocal.setZero();

        // construct the local [M] matrix
        // since [M] is symmetric, we only loop over the upper-half matrix
        for(unsigned int l = 0 ; l < elmProp.nSF*(elmProp.nSF+1)/2 ; ++l)
        {

            // sum over the GP
            double sum = 0.0;
            for(unsigned int k = 0 ; k < elmProp.nGP ; ++k)
            {
                // M_ij = sum_k{w_k*l_i(x_k)*l_j(x_k)*det[J](x_k)}
                sum += prodFunc[k][l]*mesh.elements[elm].determinantHD[k];
            }

            MLocal(IJ[l].first, IJ[l].second) = sum;

            // if we are not on the diagonal, we also add the lower-half matrix
            if(IJ[l].first != IJ[l].second)
            {
                MLocal(IJ[l].second, IJ[l].first) = sum;
            }
        }

        // inverse local M matrix (which is also symmetric)
        MLocal = MLocal.inverse();

        // set the indices for the global [M] matrix
        for(unsigned int l = 0 ; l < elmProp.nSF*(elmProp.nSF+1)/2 ; ++l)
        {
            index.push_back(Eigen::Triplet<double>
                (IJ[l].first + offsetMatrix,
                    IJ[l].second + offsetMatrix,
                    MLocal(IJ[l].first, IJ[l].second)));

            // if we are not on the diagonal, we also add the lower-half matrix
            if(IJ[l].first != IJ[l].second)
            {
                index.push_back(Eigen::Triplet<double>
                    (IJ[l].second + offsetMatrix,
                        IJ[l].first + offsetMatrix,
                        MLocal(IJ[l].first, IJ[l].second)));
            }
        }

        // increase the offset of the local matrix
        offsetMatrix += elmProp.nSF;
    }

    // add the triplets in the sparse matrix
    invM.setFromTriplets(index.begin(), index.end());
}
