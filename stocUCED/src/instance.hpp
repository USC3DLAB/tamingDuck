
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
	bool initialize (PowSys *powSys, StocProcess *stoc);
	
	PowSys		*powSys;
	StocProcess	*stoc;
	Solution	solution;
	
	ScenarioType DA_observ;
	ScenarioType ST_observ;
	ScenarioType RT_observ;

	vector< vector<double> > DA_load;	// load forecast for the DA-UC problem
	vector< vector<double> > ST_load;	// load forecast for the ST-UC problem
	vector< vector<double> > RT_load;	// real-time load for the ED problem

private:
	bool readLoadData (string filepath, vector<vector<double>> &load);
};

#endif
