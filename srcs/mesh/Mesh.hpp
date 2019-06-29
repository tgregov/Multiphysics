#ifndef Mesh2D_hpp_included
#define Mesh2D_hpp_included

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <Eigen/Sparse>

/**
 * \struct Edge
 * \brief Represents an edge.
 */
struct Edge
{
    int edgeTag;                        /**< Element tag of the edge*/

    std::vector<double> determinantLD;  /**< Determinant of the variable change for
                                            Gauss integration,
                                            evaluated at each Gauss point*/

    std::vector<int> nodeTags;          /**< Node tags for each point of the edge*/
    std::vector<std::vector<double>> nodeCoordinate; /**< Coordinate of the node*/
    std::vector<double> normal;   /**< Edge normal point outwards the element*/

    double length;

    std::pair<unsigned int, unsigned int> edgeInFront
        = std::pair<unsigned int, unsigned>(-1, -1); /**< Element index and edge
                                                        index to find the
                                                        second location of this edge
                                                         in the mesh*/
    std::vector<unsigned int> nodeIndexEdgeInFront;  /**< Track if the node are
                                                     inverted in the edge in front*/
    std::string bcName; /**< Name of the physical group of the boundary condition
                             in which the edge is in (if any)*/

    std::vector<unsigned int> offsetInElm;  /**< Offset of the edge nodes in the element*/
};

/**
 * \struct Element
 * \brief Represents an element.
 */
struct Element
{
    int elementTag;                     /**< element tag of the element*/
    int elementTypeHD;                  /**< Type of the element*/
    int elementTypeLD;                  /**< Type of the edges*/ //Store it only once

    unsigned int offsetInU;             /**< Offset of the element in the unknowns
                                            vector*/

    std::vector<double> determinantHD;  /**< Determinant of the variable change for
                                             Gauss integration, evaluated at each
                                             Gauss point*/

    std::vector<double> jacobianHD;     /**< Jacobian of the variable change for
                                            Gauss integration, evaluated at each
                                            Gauss point*/

    std::vector<Edge> edges;            /**< List of edge which compose the element */
    std::vector<int> nodeTags;          /**< List of node tags of the element */
    std::vector<std::vector<double>> nodesCoord;    /**< Node coordinates
                                                        of the elemnts */

    std::vector<Eigen::SparseMatrix<double>> dM;    /**< Partial M matrix */
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
    std::vector<std::vector<double>> prodFunc;  /**< Cross-product w_k*l_i*l_j
                                                        evaluated at each GP */

    std::vector<std::vector<double>> pondFunc;  /**< Cross-product w_k*l_i
                                                        evaluated at each GP */

    std::vector<std::pair<unsigned int, unsigned int>> IJ;  /**< Index list of
                                                                the elements of
                                                                prodFunc */
    std::vector<std::vector<double>> lalb;

};


struct NodeData
{
    unsigned int numNodes;
    std::vector<int> elementTags;
    std::vector<unsigned int> elementNumNodes;
    std::vector<int> nodeTags;
    std::vector<std::vector<double>> coord;
};


/**
 * \struct Mesh
 * \brief Represents a mesh.
 */
struct Mesh
{
    std::map<int, ElementProperty> elementProperties; /**< Store element
                                                            properties for each
                                                            element type */

    std::map<std::string, std::vector<int>> nodesTagBoundary;/**< Tags of the nodes
                                                                    per BC */

    int entityTagHD;                    /**< Tag of the HD entity*/
    int entityTagLD;                    /**< Tag of the LD entity linked to this
                                                HD entity*/

    std::vector<Element> elements;      /**< List of the elements inside the entity*/

    unsigned short dim;             /**< Mesh dimension (1, 2, (3)) */

    NodeData nodeData;

    double DxMin;
};

/**
 * \brief Get the tags of nodes (i.e. of unknowns) given a mesh.
 * \param mesh The structure that contains the mesh.
 * \return Vector containing all the tags of the nodes inside a mesh.
 */
std::vector<int> getTags(const Mesh& mesh);


/**
 * \brief Read a mesh from a file.msh
 * \param mesh The structure which will contain loaded informations.
 * \param fileName The name of the file to load.
 * \param intScheme Integration scheme for the basis functions evaluation.
 * \param basisFuncType The type of basis function you will use.
 */
bool readMesh(Mesh& mesh, const std::string& fileName,
              const std::string& intScheme, const std::string& basisFuncType);

#endif // Mesh2D_hpp_included

