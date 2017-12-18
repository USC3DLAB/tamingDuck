
/*************************************************
 
 Instance for a power production planning problem
 
 ************************************************/

#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include <iostream>
#include <string>

#include "PowSys.hpp"
#include "stoc.hpp"
#include "config.hpp"
#include "solution.hpp"

using namespace std;

class instance {

public:
    instance ();										
	void initialize (PowSys *powSys, scenarios *stoc, string inputDir);
	
	PowSys		*powSys;
	scenarios	*stoc;
	Solution	solution;
	
	vector< vector<double> > DA_load;	// load forecast for the DA-UC problem
	vector< vector<double> > ST_load;	// load forecast for the ST-UC problem
	vector< vector<double> > RT_load;	// real-time load for the ED problem
	
	// TODO: I think runParam should be here

private:
	bool readLoadData (string filepath, vector<vector<double>> &load);
};

#endif
