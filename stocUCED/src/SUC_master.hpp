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
#include "UCmodel.hpp"
#include "instance.hpp"
#include "config.hpp"
#include "misc.hpp"

class SUCmaster {

public:
	SUCmaster ();
	~SUCmaster();
	
	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep);

	bool solve ();
	double getObjValue();
	
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

	SUCsubprob	sub;		// Benders' subproblem
	UCmodel		warmUpProb;	// MIP Model for warm-starting Benders' decomposition
	void		setWarmUp();
	
	void preprocessing();

	bool getGenState(int genId, int period);				// reads from Solution.x
	void setGenState(int genId, int period, double value);	// writes to Solution.x

	vector<int>		minUpTimePeriods;	// minimum uptime in periods (obeying assumptions)
	vector<int>		minDownTimePeriods;	// minimum downtime in periods (obeying assumptions)
	
	vector<vector<double>> maxCapacity;	// expected generator capacity
	vector<double>		   sysLoad;		// aggregated load at each period
	
	int	beginMin;				// t=0 in the model corresponds to this minute in the planning horizon
	int numGen, numLine, numBus, numPeriods, numBaseTimePerPeriod, rep;
	double periodLength;
	
	char buffer[30];
};



#endif /* SUCmaster_hpp */
