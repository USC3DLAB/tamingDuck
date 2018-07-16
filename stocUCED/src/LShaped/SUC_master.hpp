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
#include <set>

#include "SUC_recourse.hpp"
#include "../UCmodel.hpp"
#include "../instance.hpp"
#include "../config.hpp"
#include "../misc.hpp"
#include "LShapedCallback.hpp"

class SUCmaster {
	friend class LShapedCallback;
public:
	SUCmaster ();
	~SUCmaster();
	
	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep);

	bool solve ();
	double getObjValue();
	
	/* Lazy-constraint callback */
	class LazySepCallbackI : public IloCplex::LazyConstraintCallbackI {
		SUCmaster & me;
	public:
		IloCplex::CallbackI* duplicateCallback() const { return (new(getEnv()) LazySepCallbackI(*this)); }
		LazySepCallbackI(IloEnv env, SUCmaster &xx) : IloCplex::LazyConstraintCallbackI(env), me(xx) {}
		void main();
	};
	
	/* Incumbent callback */
	class IncCallbackI : public IloCplex::IncumbentCallbackI {
		SUCmaster & me;
	public:
		IloCplex::CallbackI* duplicateCallback() const { return (new(getEnv()) IncCallbackI(*this)); }
		IncCallbackI(IloEnv env, SUCmaster &xx) : IloCplex::IncumbentCallbackI(env), me(xx) {}
		void main();
	};

	map< vector<bool>, vector<int> > evaluatedSolns;		// keep track of evaluated solutions, so that you don't test them again!
	instance*	inst;

private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	
	IloArray<IloNumVarArray> s, x, z, p, L, O;
	IloNumVarArray	eta;	// exp value of the 2nd-stage subproblem
	
	IloArray<IloNumArray> xvals;	// values to be passed to the subproblems
	IloRangeArray BendersCuts;		// Benders Cuts collected from the LP relaxation
	bool	LinProgRelaxFlag, MeanProbFlag;
	double	LinProgRelaxObjVal = 0;
	int		LinProgRelaxNoObjImp = 0;
	
	vector<int> rndPermutation;
	
	vector<double> expInitGen;
	
	bool multicut;
	
	
	ProblemType probType;
	ModelType	modelType;

	SUCrecourse	recourse;	// Benders' subproblem
	UCmodel		warmUpProb;	// MIP Model for warm-starting Benders' decomposition
	void		setWarmUp();
	
	void preprocessing();

	double	getEDGenProd (int genId, int period);				// reads from Solution.gED
	bool	getGenState  (int genId, int period);				// reads from Solution.x
	void	setGenState  (int genId, int period, double value);	// writes to Solution.x
	void	setUCGenProd (int genId, int period, double value);	// writes to Solution.gUC
	double	getDAUCGenProd (int genId, int period);				// reads from Solution.gUC
	double	getGenProd   (int g, int t); // reads from Solution.gED, or gUC, and handles the beginning of the planning horizon
	
	int	checkShutDownRampDownInconsistency (int g);
	
	vector<double>	minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<int>		minUpTimePeriods;	// minimum uptime in periods (obeying assumptions)
	vector<int>		minDownTimePeriods;	// minimum downtime in periods (obeying assumptions)
	
	vector<vector<double>> busLoad;		// load at each bus and period
	vector<vector<double>> expCapacity;	// expected generator capacity
	vector<double>		   sysLoad;		// aggregated load at each period
	
	int	beginMin;				// t=0 in the model corresponds to this minute in the planning horizon
	int numGen, numLine, numBus, numPeriods, numBaseTimePerPeriod, rep;
	double periodLength;
	
	char buffer[30];
};

#endif /* SUCmaster_hpp */
