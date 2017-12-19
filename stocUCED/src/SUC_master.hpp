//
//  master.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef SUCmaster_hpp
#define SUCmaster_hpp

#include <stdio.h>
#include <ilcplex/ilocplex.h>

#include "SUC_subprob.hpp"
#include "instance.hpp"
#include "config.hpp"
#include "misc.hpp"

class SUCmaster {

public:
	SUCmaster ();
	~SUCmaster();
	
	void formulate (instance &inst, int model_id);
	bool solve ();
	
//    void warmStart (Solution &soln);
//    void warmStart (vector<Solution> &soln);
	
	/* Lazy-constraint callback */
	class LazySepCallbackI : public IloCplex::LazyConstraintCallbackI {
		SUCmaster & me;
	public:
		IloCplex::CallbackI* duplicateCallback() const { return (new(getEnv()) LazySepCallbackI(*this)); }
		LazySepCallbackI(IloEnv env, SUCmaster &xx) : IloCplex::LazyConstraintCallbackI(env), me(xx) {}
		void main();
	};

private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	
	IloArray<IloNumVarArray> s, x, z;
	IloNumVar	eta;	// exp value of the 2nd-stage subproblem
	
	instance*	inst;
	
	void preprocessing();
	
//	subprob		sub;
	
	ProblemType probType;
	ModelType	modelType;
	
	int numGen, numLine, numBus, numPeriods, numBaseTimePerPeriod;
	double periodLength;
};



#endif /* master_hpp */
