#ifndef field_hpp_included
#define field_hpp_included

#include <vector>
#include <Eigen/Dense>

struct CompleteField
{
    std::vector<Eigen::VectorXd> u;

	std::vector<std::vector<Eigen::VectorXd>> flux;

	CompleteField(unsigned int numNodes, unsigned short numUnknown, unsigned short dim)
	{
        // resize each field
	    flux.resize(dim);
	    for(unsigned short i = 0 ; i < dim ; ++i)
        {
            flux[i].resize(numUnknown);
        }

		u.resize(numUnknown);

		for(unsigned short i = 0 ; i < numUnknown ; ++i)
        {
            u[i].resize(numNodes);
            for(unsigned short j = 0 ; j < dim ; ++j)
            {
                flux[j][i].resize(numNodes);
            }
        }
	}
};

/**
 * \struct Field
 * \brief Structure that contains the main unknowns & variables for the DG method.
 */
struct Field
{
	std::vector<Eigen::VectorXd> u; /**< Solution fields (of size equal to the number of scalar unknowns) */

	std::vector<std::vector<Eigen::VectorXd>> flux; /**< Physical flux fields (of size equal to the number of dimension, with each
                                                         dimension as a size equal to the number of scalar unknowns) */
    std::vector<Eigen::VectorXd> s;

	std::vector<Eigen::VectorXd> DeltaU; /**< Time-integration increment */

    std::vector<Eigen::VectorXd> Iu; /**< RHS fields */

    std::vector<Eigen::VectorXd> k1; /**< Temporary integration variables (useful for RK schemes) */
    std::vector<Eigen::VectorXd> k2; /**< Temporary integration variables (useful for RK schemes) */
    std::vector<Eigen::VectorXd> k3; /**< Temporary integration variables (useful for RK schemes) */
    std::vector<Eigen::VectorXd> k4; /**< Temporary integration variables (useful for RK schemes) */

	/**
     * \brief Constructor
     * \param numNodes The number of nodes in the mesh.
     * \param numUnknown The number of unknowns of the problem.
     * \param dim The dimension of the mesh.
     */
	Field(unsigned int numNodes, unsigned short numUnknown, unsigned short dim)
	{
        // resize each field
	    flux.resize(dim);
	    for(unsigned short i = 0 ; i < dim ; ++i)
        {
            flux[i].resize(numUnknown);
        }

		u.resize(numUnknown);
		s.resize(numUnknown);
		DeltaU.resize(numUnknown);
		Iu.resize(numUnknown);
		for(unsigned short i = 0 ; i < numUnknown ; ++i)
        {
            u[i].resize(numNodes);
            DeltaU[i].resize(numNodes);
            Iu[i].resize(numNodes);
            s[i].resize(numNodes);
            for(unsigned short j = 0 ; j < dim ; ++j)
            {
                flux[j][i].resize(numNodes);
                flux[j][i] = Eigen::VectorXd::Zero(numNodes);
            }
        }

        k1.resize(numNodes);
        k2.resize(numNodes);
        k3.resize(numNodes);
        k4.resize(numNodes);
	}
};

/**
 * \struct PartialField
 * \brief Structure that contains temporary unknowns while computing th fluxes.
 */
struct PartialField
{
    std::vector<Eigen::VectorXd> partialIu;         /**< Partial RHS */
    std::vector<std::vector<Eigen::VectorXd>> g;    /**< Partial fields (useful for the computation of the flux at the nodes) */

    std::vector<std::vector<double>> FluxAtBC;  /**< Boundary fluxes */
    std::vector<double> uAtBC;                  /**< Boundary fields (useful for the computation of the flux at BC) */

    /**
     * \brief Constructor
     * \param numUnknown The number of unknowns of the problem.
     * \param dim The dimension of the mesh.
     */
    PartialField(unsigned short numUnknown, unsigned short dim)
    {
        partialIu.resize(numUnknown);
        uAtBC.resize(numUnknown);
        g.resize(dim);
        FluxAtBC.resize(dim);
        for(unsigned short i = 0 ; i < dim ; ++i)
        {
            FluxAtBC[i].resize(numUnknown);
            g[i].resize(numUnknown);
        }
    }
};

#endif /* field_hpp */
