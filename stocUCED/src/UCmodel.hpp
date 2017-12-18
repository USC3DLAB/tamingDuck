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

public:
	UCmodel ();
	~UCmodel();

	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin);
	bool solve ();
	
private:
	/* cplex objects */
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;

	IloArray<IloNumVarArray> s, x, z, p, p_var;

	/* data */
	instance*	inst;
	ProblemType probType;
	ModelType	modelType;
	
	int	beginMin;				// t=0 in the model corresponds to this minute in the planning horizon
	int numTimePerPeriod;		// e.g., this many components in Solution object will be set by a single-period decision of the UCmodel
								// it returns, the how many ED periods are there, in a single UC period.

	int numGen;					// copied from instance->PowSys for convenience
	int numBus;					// ..
	int numLine;
	int numPeriods;
	
	double periodLength;		// in minutes

	void preprocessing ();		// see the implementation

	bool getGenState(int genId, int period);				// reads from Solution.x
	void setGenState(int genId, int period, double value);	// writes to Solution.x
	
	vector<double>	minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<int>		minUpTimePeriods;	// minimum uptime in periods (obeying assumptions)
	vector<int>		minDownTimePeriods;	// minimum downtime in periods (obeying assumptions)

	vector<vector<double>> expCapacity;	// expected generator capacity
	vector<vector<double>> busLoad;		// load at each bus and period
	vector<double>		   sysLoad;		// aggregated load at each period
	
	/* miscellaneous */
	char buffer[30];
};
#endif /* UCmodel_hpp */
