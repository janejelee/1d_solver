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

		std::cout << 1.234e5 << std::endl;
		double f,g ;

		 double dim = 2;
		 static const unsigned int temp_component           = 2+2*dim+1;
		   static const unsigned int poro_component           = 2+2*dim;

		    	std::vector<double> W(10,10);
		    	W[temp_component] = 100.0;
		    	W[poro_component] = 0.7;

		f = -9.9e-6*W[temp_component]*W[temp_component]*W[temp_component] + 0.002712*W[temp_component]*W[temp_component] -
		    		  	  	  0.8176*W[temp_component] + 1032.0;


		g =  -0.001361*W[temp_component]*W[temp_component] + 1.478*W[temp_component] +737.0;

		double kappar, kappaf, kappa;
		 kappaf = -6.0e-6*W[temp_component]*W[temp_component] + 0.00162*W[temp_component] + 0.5643;
		 kappar = 5.0e-6*W[temp_component]*W[temp_component] - 0.00585*W[temp_component] + 3.515;
		 kappa = pow(kappaf,W[poro_component]) * pow(kappar,(1.0-W[poro_component]));


		std::cout << f << std::endl << g << std::endl << kappa << std::endl;

	return 0;
}

