#ifndef Mesh2D_hpp_included
#define Mesh2D_hpp_included

#include <map>
#include <string>
#include <utility>
#include <vector>

struct Edge
{
    int edgeTag;
    std::pair<int, int> nodeTags;
};

struct Element2D
{
    int elementTag;

    std::vector<Edge> edges;
    std::vector<std::pair<double, double>> edgesNormal; //Here because outward of the element
};

struct Entity2D //problem: an entity can contain multiple elTypes ! reuse a map of jacobians, determinant...
{
    int entityTag2D;
    int entityTag1D; //Store the entity id of the edges only!

    std::vector<Element2D> elements;

    std::map<int, std::vector<double>> determinant2D;
    std::map<int, std::vector<double>> jacobian2D;

    //From the id of the Edges
    std::map<int, std::vector<double>> determinant1D;
    std::map<int, std::vector<double>> jacobian1D;
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
