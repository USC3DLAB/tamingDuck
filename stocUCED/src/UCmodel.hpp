//
//  UCmodel.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 11/11/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//
/*
#ifndef UCmodel_hpp
#define UCmodel_hpp

#include "misc.hpp"

#include <ilcplex/ilocplex.h>

#include "config.hpp"
#include "solution.hpp"
#include "instance.hpp"

class UCmodel {

public:
	UCmodel ();
	~UCmodel();

	void formulate (instance &inst, ProblemType prob_type, ModelType model_type, int first_hour);
	bool solve ();
	void exportModel ();
	void printSolution ();

	void updateSoln (solution &soln);

	solution soln;
	vector<solution> solnPool;

private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;

	IloArray <IloNumVarArray> s, x, z, p, p_var;

	instance*	inst;

	/* Model dependent data *
	void preprocessing ();

	int	nb_periods;
	int	period_len;
	int	begin_hour;

	ProblemType prob_type;

	vector< vector<double> >	demand;
	vector< vector<double> >	capacity;
	vector<double>				aggregated_demand;
	vector<double>				minimum_production_req;
	vector<int>					minimum_uptime_periods;
	vector<int>					minimum_dotime_periods;
};
#endif /* UCmodel_hpp */
