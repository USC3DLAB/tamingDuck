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

#include "../misc.hpp"

using namespace std;

struct oneStocProc {
	string	name;				// Type of random variable
	string	type;				// A descriptor
	map<string, int> mapVarNamesToIndex;	// Name of the variables
	
	vector<string> rowNames;	// Name of the rows
	int 	numVars;			// Number of random variables
	int 	numT;				// Number of time periods
	vector <vector <double> > vals;	// A two-dimensional matrix holding observation of the stochastic process
};

struct scenarioType{
	vector<string> name;				// lists all the stochastic processes for which scenarios were generated
	int numOmega;						// Number of random variables in the stochastic process
	int T;								// Time horizon
	int cnt;							// Number of scenarios
	vector<vector<vector<double>>> vals;	// scenario/observations of random variables: dim1 (rows) - time, dim2 (columns) - random variables, and dim (3) outcomes.
};

class stocProcess {

public:
	stocProcess();
	~stocProcess();
	stocProcess(string inputDir, string sysName);

	int numStocProc;
	vector <oneStocProc> sp;		// A vector of scenType with elements corresponding to each type of randomness

private:
	oneStocProc read(string fname, char delimiter, bool readColNames, bool readRowNames);
};

void read(string fname, char delimiter, bool colNames, bool rowNames);
scenarioType createScenarioList(stocProcess S, vector<int> S_indices, int T, int numVals);

#endif /* STOC_HPP_ */
