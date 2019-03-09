#ifndef Mesh2D_hpp_included
#define Mesh2D_hpp_included

#include <map>
#include <string>
#include <utility>
#include <vector>

/**
 * \struct Edge
 * \brief Represents an edge.
 */
struct Edge
{
    int edgeTag;                        /**< 1D element tag of the edge*/

    std::vector<double> determinant1D;  /**< Determinant of the variable change for Gauss integration,
                                            evaluated at each Gauss point*/

    std::pair<int, int> nodeTags;       /**< Node tags for each point of the edge*/
};

/**
 * \struct Element2D
 * \brief Represents a 2D element.
 */
struct Element2D
{
    int elementTag;                     /**< element tag of the element*/
    int elementType2D;                  /**< 2D type of the element*/
    int elementType1D;                  /**<12D type of the element*///Store it only once

    std::vector<double> determinant2D;  /**< Determinant of the variable change for Gauss integration,
                                             evaluated at each Gauss point*/
    std::vector<double> jacobian2D;     /**< Jacobian of the variable change for Gauss integration,
                                             evaluated at each Gauss point*/

    std::vector<Edge> edges;            /**< List of edge which compose the element */
    std::vector<std::pair<double, double>> edgesNormal;     /**< List of edge normal point outwards the element*/
};

/**
 * \struct Entity2D
 * \brief Represents a 2D entity.
 */
struct Entity2D
{
    int entityTag2D;                    /**< Tag of the 2D entity*/
    int entityTag1D;                    /**< Tag of the 1D entity linked to this 2D entity*/

    std::vector<Element2D> elements;    /**< List of the elements inside the entity*/

    std::map<int, std::vector<int>> elementTags2D;      /**< Tag of the element inside the entity per element type */
    std::map<int, std::vector<int>> nodesTags2D;        /**< Tag of the nodes inside the entity per element type */
    std::map<int, std::vector<int>> nodesTagsPerEdge2D; /**< Tag of the nodes per edge inside the entity per element type */
};

/**
 * \struct ElementProperty
 * \brief Store the different properties of a certain element type.
 */
struct ElementProperty
{
    std::string name;   /**< Name of the element type */
    int dim;            /**< Dimension of the element type */
    int order;          /**< Order of the element type */
    int numNodes;       /**< Number of nodes of the element type */

    std::vector<double> paramCoord;
    std::vector<double> basisFunc;      /**< Basis functions evaluated at each Gauss point in the reference axis */
    std::vector<double> basisFuncGrad;  /**< Basis functions gradient evaluated at each Gauss point in the reference axis */
    std::vector<double> intPoints;      /**< Integration points for Gauss integration */
    int numComp;
};

/**
 * \struct Mesh2D
 * \brief Represents a 2D mesh.
 */
struct Mesh2D
{
    std::map<int, ElementProperty> elementProperties2D; /**< Store element properties for each 2D element type */
    std::map<int, ElementProperty> elementProperties1D; /**< Store element properties for each 1D element type */

    std::vector<Entity2D> entities; /**< List of entities inside the mesh */
};

/**
 * \brief Read a mesh from a file.msh
 * \param mesh2D The structure which will contain loaded informations.
 * \param fileName The name of the file to load.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 * \param basisFuncGradType The type of basis function you will use ("Grad" prefix)(will be droped).
 */
bool readMesh2D(Mesh2D& mesh2D, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType,
                const std::string& basisFuncGradType);

#endif // Mesh2D_hpp_included
