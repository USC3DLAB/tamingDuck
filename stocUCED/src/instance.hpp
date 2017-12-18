
/*************************************************
 
 Instance for a power production planning problem
 
 ************************************************/

#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include <iostream>

#include "PowSys.hpp"
#include "stoc.hpp"
#include "config.hpp"
#include "solution.hpp"

using namespace std;

class instance {

public:
    instance ();										
	void initialize (PowSys *powSys, scenarios *stoc);
	
	PowSys		*powSys;
	scenarios	*stoc;
	Solution	solution;
	
	// TODO: I think runParam should be here
};

#endif
