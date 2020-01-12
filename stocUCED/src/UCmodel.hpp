//
//  UCmodel.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 11/11/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef UCmodel_hpp
#define UCmodel_hpp

#include <ilcplex/ilocplex.h>

#include "misc.hpp"
#include "config.hpp"
#include "instance.hpp"
#include "solution.hpp"

class UCmodel {
	friend class SUCmaster;		// SUCmaster uses UCmodel for warm up purposes
	
public:
	UCmodel ();
	~UCmodel();

	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep, int t0);
	double getObjValue();
	bool solve ();
	bool solve (bool saveSolution);
	
private:
	/* cplex objects */
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;

	IloArray<IloNumVarArray> s, x, z, p, p_var, L, O, v, I;
	void initializeVariables();

	/* data */
	instance*	inst;
	ProblemType probType;
	ModelType	modelType;
	
	int	beginMin;				// t=0 in the model corresponds to this minute in the planning horizon
	int numBaseTimePerPeriod;	// e.g., this many components in Solution object will be set by a single-period decision of the UCmodel
								// it returns, the how many ED periods are there, in a single UC period.

	int numGen;					// copied from instance->PowSys for convenience
	int numBus;					// ..
	int numLine;
	int numPeriods;
	int numBatteries;
	int rep;
	
	double periodLength;		// in minutes

	void preprocessing ();		// see the implementation

	bool	getGenState(int genId, int period);					// reads from Solution.x
	void	setGenState(int genId, int period, double value);	// writes to Solution.x
	void	setUCGenProd(int genId, int period, double value);	// writes to Solution.gUC
	void 	setBtState(int btId, int period, double value);		// writes to Solution.btState_UC
	double	getBatteryState (int genId, int period);			// reads from Solution.btState_ED
	double	getEDGenProd(int genId, int period);				// reads from Solution.gED
	double 	getUCGenProd(int genId, int period);				// reads from Solution.gUC

  double	getGenProd(int g, int t);		// reads from Solution.gED, or gUC, and handles the beginning of the planning horizon
	
	void 	saveSolution();
	
	int	checkShutDownRampDownInconsistency (int g);
	
	vector<double>	minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<int>		minUpTimePeriods;	// minimum uptime in periods (obeying assumptions)
	vector<int>		minDownTimePeriods;	// minimum downtime in periods (obeying assumptions)

	vector<vector<double>> capacity;	// generator capacity
	vector<vector<double>> busLoad;		// load at each bus and period
	vector<double>		   sysLoad;		// aggregated load at each period
	
	/* miscellaneous */
	char buffer[30];
};
#endif /* UCmodel_hpp */
