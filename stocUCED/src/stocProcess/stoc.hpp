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

#include <Rinside.h>

#include "../config.hpp"
#include "../misc.hpp"

using namespace std;

struct OneStocProc {
	string	name;							// Type of random variable
	string	type;							// A descriptor
	map<string, int> mapVarNamesToIndex;	// Name of variables and corresponding indices
	vector<string> rowNames;				// Name of the rows
	int 	numVars;						// Number of random variables
	int 	numT;							// Number of time periods
	vector <vector <double> > vals;			// A two-dimensional matrix holding observation of the stochastic process
};

struct ScenarioType{						// set of scenarios
	string name;							// A descriptor of the type of scenarios contained.
	map<string, int> mapVarNamesToIndex;	//
	int numVars;							// Number of random variables in the stochastic process
	int T;									// Time horizon
	int cnt;								// Number of scenarios
	vector<vector<vector<double> > > vals;	// scenario/observations of random variables: dim1 (rows) - time, dim2 (columns) - random variables, and dim (3) outcomes.
};

class StocProcess {

public:
	StocProcess();
	~StocProcess();
	StocProcess(string inputDir, string sysName);

	int numStocProc;
	vector <OneStocProc> sp;		// A vector of scenType with elements corresponding to each type of randomness	
	map<string, vector<int> > mapTypeToIndex;

private:
	OneStocProc read(string fname, char delimiter, bool readColNames, bool readRowNames);
};

void read(string fname, char delimiter, bool colNames, bool rowNames);
ScenarioType createScenarioList(StocProcess *S, vector<int> S_indices, int lenT, int stepSize, int &numVals);
ScenarioType createScenarioList(RInside &R, bool fitModel, string &dataFolder, vector<string> &stocElems, int lenT, int &numScen);

#endif /* STOC_HPP_ */
