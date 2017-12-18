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

	void formulate (instance &inst, ProblemType prob_type, ModelType model_type, int first_hour);
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
	
	int	beginMin;				// t=0 in the model corresponds to this minute in the planning horizon
	int numSolnCompsPerPeriod;	// this many components in Solution object will be set by a single-period decision of the UCmodel

	int numGen;					// copied from instance->PowSys for convenience
	int numBus;					// ..
	int numLine;
	int numPeriods;
	
	double periodLength;		// in minutes

	void preprocessing ();		// see the implementation

	bool getGenState(int genId, int period);				// reads from Solution.x
	void setGenState(int genId, int period, double value);	// writes to Solution.x
	
	vector<vector<double>> expCapacity;	// expected generator capacity
	vector<double>	minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<int>		minUpTimePeriods;	// minimum uptime in periods (obeying assumptions)
	vector<int>		minDownTimePeriods;	// minimum downtime in periods (obeying assumptions)
	
	/* miscellaneous */
	char buffer[30];
	 
	 //TODO: Clean up these
	vector< vector<double> >	demand;
	vector<double>				aggregated_demand;
};
#endif /* UCmodel_hpp */
