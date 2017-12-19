//
//  SUC_subprob.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//
/*
#ifndef subprob_hpp
#define subprob_hpp

#include <stdio.h>
#include <ilcplex/ilocplex.h>
#include <vector>

#include "instance.hpp"
#include "commons.h"

class subprob {
	friend class master;
	friend class LazySepCallbackI;
	
public:
	subprob ();
	~subprob ();
	
	void formulate (instance &inst, int model_id);
	bool solve ();
	void setMasterSoln		(vector< vector<bool> > & gen_stat);
    
    double computeLowerBound();          // a lower bound on the subproblem
    
	double getRecourseObjValue();
	
private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	
	IloRangeArray	cons;
	IloNumArray		duals;
	
	instance*	inst;
	
	int	 model_id;
	void formulate_aggregate_system ();
	void formulate_nodebased_system ();
	void formulate_production ();
	
	void setup_subproblem (int &s);
	
	double recourse_obj_val;
	
	// Benders' Cut
	double pi_b;
	vector< vector<double> > pi_T;
	
	void reset_cut_coefs ();
	void update_optimality_cut_coefs (int &s);
	void get_feasibility_cut_coefs (int &s);
	
	// First-Stage solution
	vector< vector<bool> >* gen_stat;
		
	// Variables
    IloArray< IloNumVarArray > p, L;		// production, load-shedding
	
	// Miscellaneous
	map< IloInt, double > farkasMap;
	
	double timer_t;
	
};

#endif /* subprob_hpp */
