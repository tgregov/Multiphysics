#ifndef field_hpp_included
#define field_hpp_included

#include <vector>
#include <Eigen/Dense>

struct Field
{
	// solution fields
	std::vector<Eigen::VectorXd> u; //Unknown vector, size=nbre of scalar unknown

	// physical flux fields
	std::vector<std::vector<Eigen::VectorXd>> flux; //flux fector, flux.sizesize=nbre of dim (x,y), flux[i].size()= nbre of scalar unknown

	// time-integration increment
	std::vector<Eigen::VectorXd> DeltaU;

	// rhs fields
    std::vector<Eigen::VectorXd> Iu;
    std::vector<Eigen::VectorXd> partialIu;

    std::vector<std::vector<double>> FluxAtBC;
    std::vector<double> uAtBC;
    std::vector<double> uForBC;

    std::vector<std::vector<Eigen::VectorXd>> g;

    std::vector<Eigen::VectorXd> k1;
    std::vector<Eigen::VectorXd> k2;
    std::vector<Eigen::VectorXd> k3;
    std::vector<Eigen::VectorXd> k4;

	//constructor
	Field(unsigned int numNodes, unsigned short numUnknown, unsigned short dim)
	{
	    flux.resize(dim);
	    FluxAtBC.resize(dim);
        g.resize(dim);
	    for(unsigned short i = 0 ; i < dim ; ++i)
        {
            flux[i].resize(numUnknown);
            FluxAtBC[i].resize(numUnknown);
            g[i].resize(numUnknown);
        }

		u.resize(numUnknown);
		uAtBC.resize(numUnknown);
		uForBC.resize(numUnknown);
		DeltaU.resize(numUnknown);
		Iu.resize(numUnknown);
		partialIu.resize(numUnknown);

		for(unsigned short i = 0 ; i < numUnknown ; ++i)
        {
            u[i].resize(numNodes);
            DeltaU[i].resize(numNodes);
            Iu[i].resize(numNodes);
            for(unsigned short j = 0 ; j < dim ; ++j)
                flux[j][i].resize(numUnknown);
        }

        k1.resize(numNodes);
        k2.resize(numNodes);
        k3.resize(numNodes);
        k4.resize(numNodes);
	}
};

#endif /* field_hpp */
