//
//  SUC_subprob.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef SUCsubprob_hpp
#define SUCsubprob_hpp

#include <stdio.h>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <set>

#include "instance.hpp"

using namespace std;

class SUCsubprob {
	friend class SUCmaster;
	friend class LazySepCallbackI;
	
public:
	SUCsubprob ();
	~SUCsubprob ();
	
	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin);
	bool solve ();
	void setMasterSoln		(vector< vector<bool> > & gen_stat);
	double getEDGenProd		(int genId, int period);				// reads from inst->Solution.gED
	bool getGenState		(int genId, int period);				// reads from inst->Solution.x
    
//    double computeLowerBound();          // a lower bound on the subproblem
    
	double getRecourseObjValue();
	
private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	
	IloRangeArray	cons;
	IloNumArray		duals;
	
	instance*	inst;
	
	ModelType   modelType;
	ProblemType probType;
	
	void formulate_aggregate_system ();
	void formulate_nodebased_system ();
	void formulate_production ();
	
	void setup_subproblem (int &s);
	
	double recourse_obj_val;
	
	void update_optimality_cut_coefs (int &s);
	void get_feasibility_cut_coefs (int &s);
	
	// First-Stage solution
	vector< vector<bool> >* gen_stat;
		
	// Variables
    IloArray< IloNumVarArray > p, L;		// production, load-shedding
	
	// data
	void preprocessing();
	
	int numGen;
	int numLine;
	int numBus;
	int numPeriods;
	int numBaseTimePerPeriod;
	int numScen;
	int beginMin;
	double periodLength;

	ScenarioType *scenSet;
	
	vector<double> minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<double> maxCapacity;
	
	vector<vector<double>> busLoad;		// load at each bus and period
	vector<double>		   sysLoad;		// aggregated load at each period
	
	vector<double> sceProb;				// scenario probabilities

	// Miscellaneous
	char buffer[30];
	map< IloInt, double > farkasMap;
	
	
	// Benders' Cut
	struct BendersCutCoefs {
		double pi_b;
		vector< vector<double> > pi_T;
		
		void initialize(int numRows, int numCols) {
			resize_matrix(pi_T, numRows, numCols);
		};
		
		void reset() {
			pi_b = 0;
			for (int g=0; g<pi_T.size(); g++) fill(pi_T[g].begin(), pi_T[g].end(), 0.0);
		};
	};
	BendersCutCoefs cutCoefs;
};

#endif /* SUCsubprob_hpp */
