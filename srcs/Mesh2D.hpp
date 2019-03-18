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

    std::vector<int> nodeTags;          /**< Node tags for each point of the edge*/
    std::vector<std::pair<double, double>> nodeCoordinate; /**< 2D coordinate of the node*/
    std::pair<double, double> normal;   /**< Edge normal point outwards the element*/

    //std::optional :cry:
    std::pair<unsigned int, unsigned int> edgeInFront = std::pair<unsigned int, unsigned>(-1, -1); /**<
                                                        Element index and edge index to find the
                                                        second location of this edge in the mesh*/
    std::vector<unsigned int> nodeIndexEdgeInFront;        /**< Track if the node are inverted in the edge in front*/
    std::string bcName; /**< Name of the physical group of the boundary condition
                             in which the edge is in (if any)*/

    std::vector<unsigned int> offsetInElm;
};

/**
 * \struct Element2D
 * \brief Represents a 2D element.
 */
struct Element2D
{
    int elementTag;                     /**< element tag of the element*/
    int elementType2D;                  /**< 2D type of the element*/
    int elementType1D;                  /**< 1D type of the element*/ //Store it only once

    unsigned int offsetInU;

    std::vector<double> determinant2D;  /**< Determinant of the variable change for Gauss integration,
                                             evaluated at each Gauss point*/
    std::vector<double> jacobian2D;     /**< Jacobian of the variable change for Gauss integration,
                                             evaluated at each Gauss point*/

    std::vector<Edge> edges;            /**< List of edge which compose the element */
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

    unsigned int nGP;   /**< Number of GP */
    unsigned int nSF;   /**< Number of SF */
    std::vector<std::vector<double>> prodFunc;  /**< Cross-product w_k*l_i*l_j evaluated at each GP */
    std::vector<std::vector<double>> pondFunc;  /**< Cross-product w_k*l_i evaluated at each GP */
    std::vector<std::pair<unsigned int, unsigned int>> IJ;  /**< Index list of the elements of prodFunc */

};


/**
 * \struct Mesh2D
 * \brief Represents a 2D mesh.
 */
struct Mesh2D
{
    std::map<int, ElementProperty> elementProperties2D; /**< Store element properties for each 2D element type */
    std::map<int, ElementProperty> elementProperties1D; /**< Store element properties for each 1D element type */

    std::map<std::string, std::vector<int>> nodesTagBoundary;
    //std::map<std::string, std::vector<double>> coordNodesBoundary;

    std::vector<Entity2D> entities; /**< List of entities inside the mesh */
};


/**
 * \brief Get the number of nodes (i.e. of unknowns) given a mesh.
 * \param mesh2D The structure that contains the mesh.
 */
unsigned int getNumNodes(const Mesh2D& mesh2D);


/**
 * \brief Get the tags of nodes (i.e. of unknowns) given a mesh.
 * \param mesh2D The structure that contains the mesh.
 */
std::vector<int> getTags(const Mesh2D& mesh2D);


/**
 * \brief Read a mesh from a file.msh
 * \param mesh2D The structure which will contain loaded informations.
 * \param fileName The name of the file to load.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 */
bool readMesh2D(Mesh2D& mesh2D, const std::string& fileName,
                const std::string& intScheme, const std::string& basisFuncType);

#endif // Mesh2D_hpp_included
