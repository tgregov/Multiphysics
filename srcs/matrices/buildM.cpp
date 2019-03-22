#include <iostream>
#include "buildM.hpp"


void buildM(const Mesh2D& mesh, Eigen::SparseMatrix<double>& invM)
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
            ElementProperty elmProp = mesh.elementProperties2D.at(element.elementType2D);
            std::vector<std::vector<double>> prodFunc = elmProp.prodFunc;
            std::vector<std::pair<unsigned int, unsigned int>> IJ = elmProp.IJ;

            Eigen::MatrixXd MLocal(elmProp.nSF, elmProp.nSF);
            MLocal.setZero();

            // construct the local [M] matrix for the current element
            for(unsigned int l = 0 ; l < elmProp.nSF*(elmProp.nSF+1)/2 ; ++l)
            {

                // sum over the GP
                double sum = 0;
                for(unsigned int k = 0 ; k < elmProp.nGP ; ++k)
                {
                    sum += prodFunc[k][l]*element.determinant2D[k];
                }

                MLocal(IJ[l].first, IJ[l].second) = sum;
                if(IJ[l].first != IJ[l].second)
                {
                    MLocal(IJ[l].second, IJ[l].first) = sum;
                }   

            }

            // inverse local M matrix (which is also symmetric)
            MLocal = MLocal.inverse();
            for(unsigned int l = 0 ; l < elmProp.nSF*(elmProp.nSF+1)/2 ; ++l)
            {
                index.push_back(Eigen::Triplet<double>
                    (IJ[l].first + offsetMatrix,
                        IJ[l].second + offsetMatrix,
                        MLocal(IJ[l].first, IJ[l].second)));   

                if(IJ[l].first != IJ[l].second)
                {
                    index.push_back(Eigen::Triplet<double>
                        (IJ[l].second + offsetMatrix,
                            IJ[l].first + offsetMatrix,
                            MLocal(IJ[l].first, IJ[l].second)));
                }         
            }

            offsetMatrix += elmProp.nSF;
        }
    }

    // add the triplets in the sparse matrix
    invM.setFromTriplets(index.begin(), index.end());
}
