#include <iostream>
#include <string>
#include "readMesh.h"


EIGENTYPE F(const double& t, const EIGENTYPE& v, const EIGENTYPE& M, 
			const EIGENTYPE& S, const EIGENTYPE& uBC)
{
	// EIGENTYPE flux = ...
	// EIGENTYPE F = M^(-1)*(flux - S*u)

	return F;
}


bool timeInteg(	EIGENTYPE& u, const Mesh& mesh, const string& scheme, 
				const double& h, const int& nbrTimeSteps)
{

	// comments:
	// * u is a matrix such that u[i] contains the nodal data at time i*h
	//		=> to add to the argument of the function, but I don't know what 
	//			is the type in the Eigen library
	// * we have to get M and S matrices through another function, and are 
	//		determined using the mesh variable (we should probably put the
	//		functions that determine M and S in another .cpp function, since
	//		we expect them to be quite long (maybe the same for the flux used
	//		in F()
	//		N.B.: for non-linear PDE S might be dependent on u itself, but 
	//		the scheme will be slightly changed (not much) => easily generalized
	// * we have to get the BC (the IC are contained in u[0]), denoted uBC

	if(scheme.compare("RK1")) // Runge-Kutta of order 1 (i.e. explicit Euler)
	{ 

		double t = 0;
		for(int i = 1; i < nbrTimeSteps+1; i++)
		{
			// To change to use the abilities of Eigen
			u[i] = u[i-1] + F(t, u[i-1], M, S, uBC)*h;
			t += h;
		}

	} else if(scheme.compare("RK4")){ // Runge-Kutta of order 4


		// k1, k2, k3, k4 are the same dimensions as u[i] => get the best type 
		// in Eigen
		EIGENTYPE k1, k2, k3, k4;

		for(int i = 1; i < nbrTimeSteps+1; i++)
		{
			// To change to use the abilities of Eigen
			k1 = F(t, u[i-1], M, S, uBC);
			k2 = F(t + h/2, u[i-1] + h*k1/2, M, S, uBC);
			k3 = F(t + h/2, u[i-1] + h*k2/2, M, S, uBC);
			k4 = F(t + h, u[i-1] + h*k3, M, S, uBC);

			u[i] = u[i-1] + (k1 + 2*k2 + 2*k3 + k4)*h/6;
			t += h;
		}

	} else{
		std::cerr 	<< "The integration scheme " << scheme 
					<< "is not implemented" << std::endl;
		return 1;
	}

	return 0;
}