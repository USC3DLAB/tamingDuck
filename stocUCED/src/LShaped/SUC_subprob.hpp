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

#include "BendersCutCoefs.hpp"
#include "../instance.hpp"

using namespace std;

class SUCsubprob {
	friend class SUCrecourse;
	
public:
	SUCsubprob ();
	~SUCsubprob ();
	
	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep, IloArray<IloNumArray> &masterSoln, vector<vector<double>> &expCapacity);
	
	bool solve (int mappedScen, BendersCutCoefs &cutCoefs, double &objValue, vector<double> &initGen, vector<vector<double>> &btStates);
	
	void setMasterSoln ();

private:	
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	
	IloRangeArray	cons;		// constraints which will be useful for L-shaped algorithm
	IloNumArray		duals;
	
	instance*	inst;
	
	ModelType   modelType;
	ProblemType probType;
	
	void formulate_aggregate_system ();
	void formulate_nodebased_system ();
	void formulate_production ();
	
	double 	getEDGenProd (int genId, int period);				// reads from inst->Solution.gED
	bool 	getGenState	 (int genId, int period);				// reads from inst->Solution.x
	double	getUCGenProd (int genId, int period);				// reads from Solution.gUC
	double	getGenProd   (int g, int t); 	// reads from Solution.gED, or gUC, and handles the beginning of the planning horizon
	double	getBatteryState (int batteryId, int period);
	
	int	checkShutDownRampDownInconsistency (int g);	

	void	setup_subproblem (int &s);
	double	getRandomCoef (int &s, int &t, int &loc);
	
	/* Keep record of (i) period 1 generation, (ii) battery states, later to be used by the ED model */
	void getInitGen(vector<double> &initGen);
	void getBtStates(vector<vector<double>> &btStates);

	
	// fixed master solution
	IloArray<IloNumArray> *genState;
	vector<vector<double>> *expCapacity;
	
	void compute_optimality_cut_coefs	(BendersCutCoefs &cutCoefs, int &s);
	void compute_feasibility_cut_coefs	(BendersCutCoefs &cutCoefs, int &s);
	
	// Variables
    IloArray< IloNumVarArray > p, L, x, v, I;		// production, load-shedding, state-variables (latter to be fixed by the master problem)
	
	// data
	void preprocessing();
	
	int numGen, numLine, numBus, numPeriods, numBaseTimePerPeriod, beginMin, rep, numBatteries;
	double periodLength;
	
	vector<double> minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<double> maxCapacity;
	
	vector<vector<double>> busLoad;		// load at each bus and period
	vector<double>		   sysLoad;		// aggregated load at each period
	
	// Miscellaneous
	char buffer[30];
	map< IloInt, double > farkasMap;
	
	double solve_t, setup_t, cut_t;
};

#endif /* SUCsubprob_hpp */
