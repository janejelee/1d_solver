/*
 * test_numbers.cpp
 *
 *  Created on: 15 Dec 2016
 *      Author: leej
 */



#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <memory>

int main(int argc, char*argv[])
{
	double compute_bulk_density (std::vector<double>& W, double& rho_f);

		std::cout << 1.234e5 << std::endl;
		double rho_f,c_r ;

		 double dim = 2;
		 static const unsigned int temp_component           = 2+2*dim+1;
		   static const unsigned int poro_component           = 2+2*dim;

		    	std::vector<double> W(10,10);
		    	W[temp_component] = 100.0;
		    	W[poro_component] = 0.7;

		rho_f = -9.9e-6*W[temp_component]*W[temp_component]*W[temp_component] + 0.002712*W[temp_component]*W[temp_component] -
		    		  	  	  0.8176*W[temp_component] + 1032.0;
		c_r =  -0.001361*W[temp_component]*W[temp_component] + 1.478*W[temp_component] +737.0;


		double kappar, kappaf, kappa;
		 kappaf = -6.0e-6*W[temp_component]*W[temp_component] + 0.00162*W[temp_component] + 0.5643;
		 kappar = 5.0e-6*W[temp_component]*W[temp_component] - 0.00585*W[temp_component] + 3.515;
		 kappa = pow(kappaf,W[poro_component]) * pow(kappar,(1.0-W[poro_component]));


		 double mu_0 = 2.0, E = 10.0, R = 5.0, T_0 =10.0;
		 double mu_r;
		 mu_r = mu_0*(exp((E/R) * (1.0/ W[temp_component] - 1.0/T_0)));
		 // std::cout << rho_f << std::endl << c_r << std::endl << kappa << std::endl;
		 	 std::cout << mu_r <<std::endl;

		 double rho_b;
		 rho_b = compute_bulk_density(W, rho_f);

			std::cout <<  rho_f << std::endl << rho_b << std::endl;



	return 0;
}



	double compute_bulk_density (std::vector<double>& W, double& rho_f)
	        	{
		 double dim = 2;
				 static const unsigned int temp_component           = 2+2*dim+1;
				   static const unsigned int poro_component           = 2+2*dim;
				   double rho_r = 1000.0;

	        		return ( (1.0-W[poro_component])*rho_r + W[poro_component]*rho_f );

	        	}
