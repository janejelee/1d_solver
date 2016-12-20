/*
 * ex1-3.cpp
 *
 *  Created on: 11 Feb 2016
 *      Author: leej
 */

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <memory>

int main(int argc, char*argv[])
{

	      int dim =3;
	      std::vector<std::string> names (1, "rock_pressure");
	      names.push_back ("fluid_pressure");
	      std::vector<std::string> vr (dim, "rock_velocity");
	      std::vector<std::string> vf (dim, "fluid_velocity");
	      names.insert(names.end(), vf.begin(), vf.end());
	      names.insert(names.end(), vr.begin(), vr.end());
	      names.push_back ("porosity");
	      names.push_back ("temperature");
	      std::cout << names[0] << "\n" << names[1] << "\n"<<  names[2] << "\n" << names[3] << "\n" << names[4] << "\n" <<names[5] << "\n" << names[6] << "\n"<< names[7] << "\n";


	return 0;
}



