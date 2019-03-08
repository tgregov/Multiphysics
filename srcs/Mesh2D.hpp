#ifndef Mesh2D_hpp_included
#define Mesh2D_hpp_included

#include <map>
#include <string>
#include <utility>
#include <vector>

struct Edge
{
    int edgeTag;

    std::vector<double> determinant1D;

    std::pair<int, int> nodeTags;
};

struct Element2D
{
    int elementTag;
    int elementType2D;
    int elementType1D; //Store it only once

    std::vector<double> determinant2D;
    std::vector<double> jacobian2D;

    std::vector<Edge> edges;
    std::vector<std::pair<double, double>> edgesNormal; //Here because outward of the element
};

struct Entity2D //problem: an entity can contain multiple elTypes ! reuse a map of jacobians, determinant...
{
    int entityTag2D;
    int entityTag1D; //Store the entity id of the edges only!

    std::vector<Element2D> elements;

    std::map<int, std::vector<int>> elementTags2D;
    std::map<int, std::vector<int>> nodesTags2D;
    std::map<int, std::vector<int>> nodesTagsPerEdge2D;
};

struct ElementProperty
{
    std::string name;
    int dim;
    int order;
    int numNodes;

    std::vector<double> paramCoord;
    std::vector<double> basisFunc;
    std::vector<double> basisFuncGrad;
    std::vector<double> intPoints;
    int numComp;
};

struct Mesh2D
{
    std::map<int, ElementProperty> elementProperties2D;
    std::map<int, ElementProperty> elementProperties1D;

    std::vector<Entity2D> entities;
};

bool readMesh2D(Mesh2D& mesh2D, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType,
                const std::string& basisFuncGradType);

#endif // Mesh2D_hpp_included
