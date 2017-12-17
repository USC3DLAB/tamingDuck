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
/*	void exportModel ();
	void printSolution ();

	void updateSoln (solution &soln);

//	Solution soln;
//	vector<solution> solnPool;
*/
private:
	/* cplex objects */
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;

	IloArray<IloNumVarArray> s, x, z, p, p_var;

	/* data */
	instance*	inst;
	ProblemType probType;
	
	bool getGenState(int genId, int period);
	
	/* model-dependent components */
	void preprocessing ();
	
	vector<vector<double>> expCapacity;	// expected generator capacity
	vector<double>	minGenerationReq;	// minimum production requirements
	vector<int>		minUpTimePeriods;	// minimum uptime in periods
	vector<int>		minDownTimePeriods;	// minimum downtime in periods
	
	/* for convenience */
	int numGen;
	int numBus;
	int numLine;
	int numPeriods;

	double periodLength; // minutes
	
	/* miscellaneous */
	char buffer[30];
	 
	 //TODO: Clean up these
//	int	begin_hour;

	

	
	vector< vector<double> >	demand;
	vector<double>				aggregated_demand;
};
#endif /* UCmodel_hpp */
