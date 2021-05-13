//
//  LShapedCallback.hpp
//  tamingDuck
//
//  Created by Semih Atakan on 3/8/18.
//  Copyright © 2018 Semih Atakan. All rights reserved.
//

#ifndef LShapedCallback_hpp
#define LShapedCallback_hpp

#include <map>
#include <vector>
#include <stdio.h>
#include <ilcplex/ilocplex.h>
#include "../config.hpp"
#include "ilconcert/ilothread.h"		// for mutex (must be included below boost libraries to avoid errors)

using namespace std;

class SUCmaster;

class LShapedCallback: public IloCplex::Callback::Function {
public:
	// constructor
	LShapedCallback (SUCmaster &master);
	
	// destructor
	virtual ~LShapedCallback () {}
	
	// when to invoke the callback? this function will tell you
	virtual void invoke (const IloCplex::Callback::Context &context);
	
	// L-shaped cut generator
	inline bool addLShapedCuts (const IloCplex::Callback::Context &context);
	
	// set expected generation amounts
	inline void saveSubprobSolns (const IloCplex::Callback::Context &context);
	
//	// randomized-rounding heuristic
//	inline void randomizedRounding (const IloCplex::Callback::Context &context);
	
	// searches solutions in the solution list
	inline map<vector<bool>, double>::iterator findSolution(const IloCplex::Callback::Context &context);
	
private:
	// empty constructor
	LShapedCallback () {}
	
	// copy constructor
	LShapedCallback (const LShapedCallback &tocopy);
	
	// pointer to the master problem
	SUCmaster *master;
	
	// parallel thread management
	IloFastMutex mutex;
	
	// solution management
	map<vector<bool>, double> solMap;
};

#endif /* LShapedCallback_hpp */
