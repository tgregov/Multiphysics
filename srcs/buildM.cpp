#include <iostream>
#include <Eigen/Sparse>
#include "buildM.hpp"
#include "Mesh2D.hpp"


// note to myself: wtf with const& mesh and [] for elementProperties2D ?
void buildM(Mesh2D& mesh, Eigen::SparseMatrix<double>& M)
{

    // the index is a vector of triplets that contains the coordinates in the [M]
    // matrix of each ot its components; the offsetMatrix is the upper-left
    // coordinate at which the current element matrix should be added
    std::vector<Eigen::Triplet<double>> index;
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
            ElementProperty elmProp = mesh.elementProperties2D[element.elementType2D];
            std::vector<std::vector<double>> prodFunc = elmProp.prodFunc;
            std::vector<std::pair<unsigned int, unsigned int>> IJ = elmProp.IJ;

            // construct the local [M] matrix for the current element
            for(unsigned int l = 0 ; l < elmProp.nSF*(elmProp.nSF+1)/2 ; ++l)
            {

                // sum over the GP
                double sum = 0;
                for(unsigned int k = 0 ; k < elmProp.nGP ; ++k)
                {
                    sum += prodFunc[k][l]*element.determinant2D[k];
                }

                // add the sum to the index vector
                index.push_back(Eigen::Triplet<double>
                    (IJ[l].first + offsetMatrix,
                    IJ[l].second + offsetMatrix,
                    sum));

                // if we are not on the diagonal, we add the symmetric part of the
                // matrix (indeed, since [M] is symmetric, the products w_k*l_i*l_j
                // were only computed once
                if(IJ[l].first != IJ[l].second)
                {
                    index.push_back(Eigen::Triplet<double>
                        (IJ[l].second + offsetMatrix,
                          IJ[l].first + offsetMatrix,
                         sum));
                }
            }

            offsetMatrix += elmProp.nSF;
        }
    }

    // add the triplets in the sparse matrix
    M.setFromTriplets(index.begin(), index.end());
}
