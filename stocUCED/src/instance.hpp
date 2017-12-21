
/*************************************************
 
 Instance for a power production planning problem
 
 ************************************************/

#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include "misc.hpp"
#include "config.hpp"
#include "./stocProcess/stoc.hpp"
#include "./powerSys/PowSys.hpp"


#include "solution.hpp"

using namespace std;

class instance {

public:
    instance ();										
	bool initialize (PowSys *powSys, StocProcess *stoc, vector<string> stocElems, vector<string> detElems);
	bool printSolution(string filepath);
	
	
	PowSys		*powSys;
	StocProcess	*stoc;
	Solution	solution;

	vector<string> hierarchy;		// Hierarchy of information/data flow. Example, "DA", "ST", "RT".
	
	vector<string> detElems;		// list of deterministic elements in the model
	vector<string> stocElems;		// list of stochastic elements in the model

	vector<ScenarioType> detObserv;  // Deterministic observations.
	vector<ScenarioType> stocObserv; // Stochastic observations.

private:
	void summary();
};

#endif
