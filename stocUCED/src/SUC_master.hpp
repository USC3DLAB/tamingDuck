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
	
	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin);

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
	
	ProblemType probType;
	ModelType	modelType;

//	subprob		sub;
	
	void preprocessing();

	bool getGenState(int genId, int period);				// reads from Solution.x
	void setGenState(int genId, int period, double value);	// writes to Solution.x

	vector<double>	minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<int>		minUpTimePeriods;	// minimum uptime in periods (obeying assumptions)
	vector<int>		minDownTimePeriods;	// minimum downtime in periods (obeying assumptions)
	
	vector<vector<double>> expCapacity;	// expected generator capacity
	vector<vector<double>> busLoad;		// load at each bus and period
	vector<double>		   sysLoad;		// aggregated load at each period
	
	int	beginMin;				// t=0 in the model corresponds to this minute in the planning horizon
	int numGen, numLine, numBus, numPeriods, numBaseTimePerPeriod;
	double periodLength;
};



#endif /* master_hpp */
