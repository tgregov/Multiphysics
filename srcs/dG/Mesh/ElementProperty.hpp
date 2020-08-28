#pragma once
#ifndef ELEMENTPROPERTY_HPP_INCLUDED
#define ELEMENTPROPERTY_HPP_INCLUDED

namespace dG
{
    struct ElementProperty
    {
        int type                            = -1;
        std::string name                    = {};
        int order                           = -1;
        int numNodes                        = -1;
        std::vector<double> localNodeCoord  = {};
        int numPrimaryNodes                 = -1;

        std::vector<double> intPointsCoord  = {};
        std::vector<double> intPointsWeigth = {};

        int basisFuncNumComp                = -1;
        std::vector<double> basisFunctions  = {};

        int basisFuncNumCompGrad                = -1;
        std::vector<double> basisFunctionsGrad  = {};
    };
}

#endif // ELEMENTPROPERTY_HPP_INCLUDED
