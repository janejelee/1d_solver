/*
 * test_switch.cpp
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

		 const int dim = 4;
		    static const unsigned int n_components             = 2*dim + 4;
		    static const unsigned int VES_component      	   = 0;
		    static const unsigned int pf_component             = 1;
		    static const unsigned int vf_first_component       = 2;
		    static const unsigned int vr_first_component       = 2+dim;
		    static const unsigned int poro_component           = 2+2*dim;
		    static const unsigned int temp_component           = 2+2*dim+1;

		    std::vector<double> forcing(n_components,1234.0);


		for (unsigned int c=0; c<n_components; c++)
			switch (c)
			{
			case vf_first_component+dim-1:
				forcing[c] = 100.0;
				break;
			case poro_component:
				forcing[c] = 0.7;
				break;
			default:
				forcing[c] = 0.0;

			}

		switch (dim)
		{
		case 2:
			forcing[vf_first_component] = 50.0;
			break;
		case 3:
			forcing[vf_first_component] = 50.0;
			forcing[vf_first_component+1] = 50.0;
			break;
		}




		 for (unsigned int n=0; n<n_components; n++)
		 {
			 std::cout << forcing[n] << std::endl;
		 }

		    	return 0;
}
