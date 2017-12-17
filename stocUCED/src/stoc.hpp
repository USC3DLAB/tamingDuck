/*
 * stoc.hpp
 *
 *  Created on: Dec 14, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef STOC_HPP_
#define STOC_HPP_

/* This header file contains all classes and structures pertaining to the stochastic processes of interest */

#include "misc.hpp"

using namespace std;

struct scenType {
	string	name;		// Type of random variable
	string	type;		// A descriptor
	vector<string> varNames;
	int 	numL;		// Number of random variables
	int 	numT;		// Number of time periods
	int 	cnt;		// Number of scenarios
	vector <vector <vector <double> > > vals;	// A three-dimensional matrix with scenario values
};

class scenarios {

public:
	scenarios();
	~scenarios();
	scenarios(string inputDir, string sysName);
	bool read(string fname, char delimiter);

	vector <scenType> scen;		// A vector of scenType with elements corresponding to each type of randomness
};


#endif /* STOC_HPP_ */
