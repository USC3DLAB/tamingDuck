
/*************************************************
 
 Instance for a power production planning problem
 
 ************************************************/

#ifndef _INSTANCE_H_
#define _INSTANCE_H_

//#include <Rinside.h>

#include "misc.hpp"
#include "config.hpp"
#include "./stocProcess/stoc.hpp"
#include "./powerSys/PowSys.hpp"
#include "solution.hpp"

using namespace std;

class instance {

public:
    instance ();
	bool initialize (PowSys *powSys, StocProcess *stoc, string RScriptsPath);
	bool printSolution(string filepath);
	
	PowSys		*powSys;
	StocProcess	*stoc;
	Solution	solution;
	
	vector<string> elems;			// list of all processes
	vector<string> detElems;		// list of deterministic elements in the model
	vector<string> stocElems;		// list of stochastic elements in the model
	
	void simulateScenarios();
	ScenarioType simulations;
	
	map< string, ScenarioType > observations;
	
	vector<ScenarioType> detObserv;  // Deterministic observations.
	vector<ScenarioType> stocObserv; // Stochastic observations.

private:
	void summary();
	
	RInside R;
};

#endif
