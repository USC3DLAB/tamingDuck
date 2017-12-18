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
	string	name;				// Type of random variable
	string	type;				// A descriptor
	vector<string> varNames;	// Name of the variables
	vector<string> rowNames;	// Name of the rows
	int 	numVars;			// Number of random variables
	int 	numT;				// Number of time periods
	vector <vector <double> > vals;	// A two-dimensional matrix with scenario values
};

struct observType{
	vector<string> name;
	int numOmega;
	int T;
	int cnt;
	vector<vector<vector<double>>> vals;
};

class scenarios {

public:
	scenarios();
	~scenarios();
	scenarios(string inputDir, string sysName);

	int numStocProc;
	vector <scenType> stocProc;		// A vector of scenType with elements corresponding to each type of randomness

private:
	scenType read(string fname, char delimiter, bool readColNames, bool readRowNames);
};

void read(string fname, char delimiter, bool colNames, bool rowNames);
observType createObservList(scenarios S, vector<int> S_indices, int T, int numVals);

#endif /* STOC_HPP_ */
