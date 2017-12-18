
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
#include "Load.hpp"

using namespace std;

class instance {

public:
    instance ();										
	void initialize (PowSys *powSys, scenarios *stoc, string inputDir);
	
	PowSys		*powSys;
	scenarios	*stoc;
	Solution	solution;
	
	Load		DA_load;	// load forecast for the DA-UC problem
	Load		ST_load;	// load forecast for the ST-UC problem
	Load		RT_load;	// real-time load for the ED problem
	
	// TODO: I think runParam should be here
};

#endif
