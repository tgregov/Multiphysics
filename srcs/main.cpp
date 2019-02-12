#include <iostream>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif // M_PI
#include <gmsh.h>

int main()
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::open("triangle.geo");
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);

    for(unsigned int i = 0 ; i < entities.size() ; ++i)
    {
        std::vector<int> nodeTags;
        std::vector<double> nodeCoords, nodeParams;
        int dim = entities[i].first, tag = entities[i].second;
        gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

        std::vector<int> elemTypes;
        std::vector<std::vector<int>> elemTags, elemNodeTags;
        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

        std::string type;
        gmsh::model::getType(dim, tag, type);
        std::cout<<"Entitites "<<i<<", type"<<type<<std::endl;
    }

    gmsh::finalize();
    return 0;
}
