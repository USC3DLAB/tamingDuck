//
//  SUC_master.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "SUC_master.hpp"

extern runType runParam;

ILOHEURISTICCALLBACK2(rounding, IloArray<IloNumVarArray> &, x, SUCmaster &, me) {
	
	/* Heuristic frequency: Every 100 nodes, after node 10 */
	if (getNnodes() < 10 || getNnodes() % 100 != 0){
		return;
	}
	
	/* Construct a Solution */
	IloNumVarArray	vars (getEnv());
	IloNumArray		vals (getEnv());
	vector<bool> 	vals_bool;
	
	for (int g=0; g<x.getSize(); g++) {
		vars.add( x[g] );
	
		IloNumArray temp (getEnv());
		getValues(temp, x[g]);
		for (int t=0; t<temp.getSize(); t++) {
			//temp[t] = round(temp[t]);
			temp[t] = (double(rand())/double(RAND_MAX) < temp[t])*1.0;
			vals_bool.push_back( round(temp[t]) );
		}
		vals.add(temp);
		temp.end();
	}
	
	/* Evaluate the Solution */
	try {
		// Important: Check if the produced solution is evaluated before. If it is, don't evaluate!
		// It is important that you don't visit the same solution multiple times, especially
		// if you're not evaluating them in the lazy-constraint callback
		auto it = me.evaluatedSolns.find(vals_bool);
		if (it == me.evaluatedSolns.end()) {
			me.inst->out() << "Heuristic is evaluating a new solution" << endl;
			setSolution(vars, vals);
			solve();
		} else {
			me.inst->out() << "Heuristic has ignored Sol" << me.evaluatedSolns[vals_bool][2] << endl;
		}
	}
	catch (IloException &e) {
		cout << e << endl;
	}

	vals.end();
	vars.end();
}

/****************************************************************************
 * Incumbent Callback
 * Used for retrieving expected initial-generation amounts for the optimal
 * master solution.
 ****************************************************************************/
IloCplex::Callback IncCallback(IloEnv env, SUCmaster &me) {
	return (IloCplex::Callback(new(env) SUCmaster::IncCallbackI(env, me)));
}

void SUCmaster::IncCallbackI::main() {
	// if the incumbent's objective is improving with the new solution
	if ( getObjValue() < getIncumbentObjValue() ) {
		
		// request expected initial generation amounts
		vector<double> expInitGen = me.recourse.getExpInitGen();
		
		// set the generation amounts
		for (int g=0; g<me.numGen; g++) {
			me.setUCGenProd(g, 0, expInitGen[g]);
		}
	}
}

/****************************************************************************
 * Lazy Constraint Callback
 * Used for managing the Benders' (L-shaped) algorithm.
 ****************************************************************************/
IloCplex::Callback LazySepCallback(IloEnv env, SUCmaster & me) {
	return (IloCplex::Callback(new(env) SUCmaster::LazySepCallbackI(env, me)));
}

bool warmup = true;
void SUCmaster::LazySepCallbackI::main()
{
	if (warmup) {
		if ( (getNodeId())._id > 0) {
			warmup = false;
			me.evaluatedSolns.clear();
		}
	}
	
	if ( me.LinProgRelaxFlag ) {
		me.inst->out() << "LinProgRelax = " << getObjValue() << endl;
		
		if ( fabs(getObjValue() - me.LinProgRelaxObjVal) / (fabs(me.LinProgRelaxObjVal)+1e-14) < 0.005 ) {
			me.LinProgRelaxNoObjImp++;
		} else {
			me.LinProgRelaxNoObjImp = 0;
		}
		
		if ( me.LinProgRelaxNoObjImp > 5 ) {
			return;
		}
		me.LinProgRelaxObjVal = getObjValue();
	}
	
	// get the solution
	for (int g=0; g<me.numGen; g++) {
		if (me.LinProgRelaxFlag) {
			getValues(me.xvals[g], me.x[g]);
		} else {
			for (int t=0; t<me.numPeriods; t++) {
				me.xvals[g][t] = round(getValue( me.x[g][t] ));
			}
		}
	}
	
	// add solution to the pool, test if it has been seen before
	vector<bool> x;
	if (!me.LinProgRelaxFlag) {
		for (int g=0; g<me.numGen; g++) {
			for (int t=0; t<me.numPeriods; t++) {
				x.push_back( round(me.xvals[g][t]) );
			}
		}
		// test if the solution was in the pool, no cut is necessary if it has already been solved before
		auto it = me.evaluatedSolns.find( x );
		if (it == me.evaluatedSolns.end()) {
			vector<int> val; val.push_back(1); val.push_back( me.evaluatedSolns.size()+1 );
			me.evaluatedSolns.insert( pair<vector<bool>, vector<int>> (x, val) );
			me.MeanProbFlag = true;

//			me.inst->out() << "(Sol " << (me.evaluatedSolns[x])[1] << ": solving mean-scenario. obj= " << getObjValue() << " eta= " << getValue(me.eta[0]) << ")" << endl;
		}
		else {
			if ( (me.evaluatedSolns[x])[0] == 1) {
				me.MeanProbFlag = false;
				(me.evaluatedSolns[x])[0] = 2;
//				me.inst->out() << "(Sol " << (me.evaluatedSolns[x])[1] << ": solving all scenarios. obj= " << getObjValue() << " eta= " << getValue(me.eta[0]) << ")" << endl;
			}
			else {
//				cout << "Returning since I've seen this soln?" << endl;
//				me.inst->out() << "(Sol " << (me.evaluatedSolns[x])[1] << ": no cut! I've seen this soln before. obj= " << getObjValue() << " eta= " << getValue(me.eta[0]) << ")" << endl;
				return;
			}
		}
	}
	
	me.evaluatedSolns[x][0] = 2;
	me.MeanProbFlag = false;
	
	// set the master solution in the subproblems
	me.recourse.setMasterSoln();
	
//	me.recourse.solveMeanProb();
//	double mean = me.recourse.getScenObjValue(0);
//	me.recourse.solve();
//	double scen = me.recourse.getObjValue();
//	if ( mean > scen ) {
//		cout << mean << "," << scen << endl;
//	}
	// solve the subproblems
	double time_t = get_wall_time();
	bool isFeasible;
	if (me.MeanProbFlag) {
		isFeasible = me.recourse.solveMeanProb();
	}
	else {
		isFeasible = me.recourse.solve();
	}
	me.inst->out() << "- Sol " << me.evaluatedSolns[x][1] << ": Obj= " << getObjValue() << " eta= " << getValue(me.eta[0]);
	if (me.MeanProbFlag) 	me.inst->out() << " RecObj= " << me.recourse.getScenObjValue(0) << " " << flush;
	else					me.inst->out() << " RecObj= " << me.recourse.getObjValue() << " " << flush;

	cout << get_wall_time() - time_t << endl;
	
	/*
	cout << me.recourse.getObjValue() << endl;
	for (int g=0; g<me.numGen; g++) {
		cout << g+1 << " " << me.xvals[g] << endl;
	}
	 */
	/* adjust for the constants
	double constant = 0;
	for (int b=0; b<me.numBus; b++) {
		constant += getValue( IloSum(me.L[b]) ) * loadShedPenaltyCoef;
		constant += getValue( IloSum(me.O[b]) ) * overGenPenaltyCoef;
	}
	for (int g=0; g<me.numGen; g++) {
		Generator *genPtr = &(me.inst->powSys->generators[g]);
		constant += getValue( IloSum(me.p[g]) ) * genPtr->variableCost*me.periodLength/60.0;
	}
	for (int s=0; s<runParam.numLSScen; s++) {
		me.recourse.objValues[s] -= constant;
	}
	 */
	
	const double AbsOptTol = 1e-6;
	const double RelOptTol = 1e-4;
	
    // check if we need to add a cut (feasible and not within optimality tolerances, or infeasible)
    bool addCut = true;
    if (isFeasible) {
		if (!me.MeanProbFlag) {
			double diff = me.recourse.getObjValue() - getValue(IloSum(me.eta));
			if ( diff/(fabs(me.recourse.getObjValue())+1e-14) < RelOptTol || diff < AbsOptTol ) {
				addCut = false;
			}
			me.inst->out() << " RecRelGap= " << diff/(fabs(me.recourse.getObjValue())+1e-14) << " RecAbsGap= " << diff << " " << endl;
		} else {
			double diff = me.recourse.objValues[0] - getValue(IloSum(me.eta));
			if ( diff/(fabs(me.recourse.getObjValue())+1e-14) < RelOptTol || diff < AbsOptTol ) {
				addCut = false;
			}
			me.inst->out() << " RecRelGap= " << diff/(fabs(me.recourse.getScenObjValue(0))+1e-14) << " RecAbsGap= " << diff << " " << endl;
		}
    }
	
	if (me.MeanProbFlag && !addCut && isFeasible) {
		me.inst->out() << "  (useless mean-scen cut) " << endl;
		//me.inst->out() << "(Sol " << (me.evaluatedSolns[x])[1] << ": solving all scenarios, mean-cut was useless. obj= " << getObjValue() << " eta= " << getValue(me.eta[0]) << ")" << endl;
		me.MeanProbFlag = false;
		me.recourse.solve();

		(me.evaluatedSolns[x])[0] = 2;
		addCut = true;
		double diff = me.recourse.getObjValue() - getValue(IloSum(me.eta));
		if ( diff/(fabs(me.recourse.getObjValue())+1e-14) < RelOptTol || diff < AbsOptTol ) {
			addCut = false;
		}
	}
	
	if (!addCut) {
		return;
	}
	
	if (isFeasible) {
		/* optimality cut */
		
		if ( me.multicut ) {
			cout << "No multicut" << endl;
//			int optcut_cnt = 0;
//			for (int s=0; s<me.eta.getSize(); s++) {
//
//				double diff = me.recourse.getScenObjValue(s) - getValue(me.eta[s]) - EPSzero;
//
//				if ( diff/(fabs(me.recourse.getScenObjValue(s))+1e-14) < me.cplex.getParam(IloCplex::EpGap)
//					|| diff < me.cplex.getParam(IloCplex::EpAGap) ) {
//					continue;
//				}
//
//				optcut_cnt++;
//
//				// get the Benders's cut coefs
//				BendersCutCoefs *cutCoefs = &(me.recourse.cutCoefs[s]);
//
//				// create the cut
//				IloExpr pi_Tx (me.env);
//				for (int g=0; g<me.numGen; g++) {
//					for (int t=0; t<me.numPeriods; t++) {
//						if ( fabs(cutCoefs->pi_T[g][t]) > 1e-10 ) {
//							pi_Tx += cutCoefs->pi_T[g][t] * me.x[g][t];
//						}
//					}
//				}
//
//				IloRange BendersCut;
//				BendersCut = IloRange(me.env, cutCoefs->pi_b, me.eta[s]-pi_Tx);	// optimality cut
//
//				// add the cut
//				try {
//					if (me.LinProgRelaxFlag) {
//						me.BendersCuts.add(BendersCut);
//						add(BendersCut);
//					}
//					else {
//						add(BendersCut).end();
//					}
//				}
//				catch (IloException &e) {
//					cout << "Exception: " << e << endl;
//					BendersCut.end();
//				}
//				pi_Tx.end();
//			}
//
//			me.inst->out() << "(optcut " << optcut_cnt << "/" << me.eta.getSize() << ")" << endl;
		}
		else
		{
			double sceProb = 1.0/(double)runParam.numLSScen;
			
			IloExpr pi_Tx (me.env);
			IloNum  pi_b = 0;
			
			for (int s=0; s<runParam.numLSScen; s++) {
				// get the Benders's cut coefs
				BendersCutCoefs *cutCoefs = &(me.recourse.cutCoefs[s]);
				
				// create the cut
				for (int g=0; g<me.numGen; g++) {
					for (int t=0; t<me.numPeriods; t++) {
						if ( fabs(cutCoefs->pi_T[g][t]) > 1e-10 ) {
							pi_Tx += cutCoefs->pi_T[g][t] * me.x[g][t];
						}
					}
				}
				pi_b += cutCoefs->pi_b;

				if (me.MeanProbFlag) break;
			}

			if (!me.MeanProbFlag) {
				pi_Tx *= sceProb;
				pi_b  *= sceProb;
			}
			
			IloRange BendersCut (me.env, pi_b, me.eta[0]-pi_Tx);	// optimality cut

			// add the cut
			try {
//				if (me.LinProgRelaxFlag || me.MeanProbFlag) {
//					me.BendersCuts.add(BendersCut);
//					add(BendersCut);
//				}
//				else {
					add(BendersCut).end();
//				}
			}
			catch (IloException &e) {
				cout << "Exception: " << e << endl;
				BendersCut.end();
			}
			pi_Tx.end();
			
			me.inst->out() << "  (optcut)" << endl;
		}
	}
	else {
		/* feasibility cut */
		me.inst->out() << "(feascut)" << endl;

		// get the Benders's cut coefs
		BendersCutCoefs *cutCoefs = &(me.recourse.cutCoefs[ me.recourse.getInfeasScenIndex() ]);
		
		// create the cut
		IloExpr pi_Tx (me.env);
		for (int g=0; g<me.numGen; g++) {
			for (int t=0; t<me.numPeriods; t++) {
				if ( fabs(cutCoefs->pi_T[g][t]) > 1e-10 ) {
					pi_Tx += cutCoefs->pi_T[g][t] * me.x[g][t];
				}
			}
		}
		
		IloRange BendersCut;
		BendersCut = IloRange(me.env, cutCoefs->pi_b, -pi_Tx);			// feasibility cut
		
		// add the cut
		try {
			add(BendersCut).end();
		}
		catch (IloException &e) {
			cout << "Exception: " << e << endl;
			BendersCut.end();
		}
		pi_Tx.end();
	}
}

SUCmaster::SUCmaster () {
	model = IloModel(env);
	cplex = IloCplex(env);
}

SUCmaster::~SUCmaster() {
	env.end();
}

/****************************************************************************
 * preprocessing
 * - Initializes basic parameters
 * - Fills the following containers, according to model units & assumptions:
 *   minUpTimePeriods
 *   minDownTimePeriods
 *   expCapacity
 *   aggLoads
 * - If it is a DA-UC problem, sets the capacity and minGenerationReq to 0
 * for generators which must be scheduled at ST-UC.
 *
 * - Assumption 1: Min generation requirement must be less than ramping
 * rates, otherwise, generators cannot switch on/off.
 * - Assumption 2: Min up/down times must at least be 1 period.
 ****************************************************************************/
void SUCmaster::preprocessing ()
{
	/* basic parameters */
	numGen	   = inst->powSys->numGen;
	numLine    = inst->powSys->numLine;
	numBus     = inst->powSys->numBus;
	
	numPeriods	 = (probType == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
	periodLength = (probType == DayAhead) ? runParam.DA_resolution : runParam.ST_resolution;
	
	numBaseTimePerPeriod = (int)round(periodLength) / runParam.ED_resolution;
	
	
	/* initialize containers */
	minGenerationReq.resize(numGen);					// minimum production requirements
	minUpTimePeriods.resize(numGen);					// minimum uptime in periods
	minDownTimePeriods.resize(numGen);					// minimum downtime in periods
	resize_matrix(expCapacity, numGen, numPeriods);		// max generator capacities
	sysLoad.resize(numPeriods);							// aggregated system load
	resize_matrix(busLoad, numBus, numPeriods);
	xvals = IloArray<IloNumArray> (env, numGen);
	for (int g=0; g<numGen; g++) {
		xvals[g] = IloNumArray (env, numPeriods);
	}
	
	
	/* Min Generation Amounts */
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		// do not remove this assignment, as STUC generators' minGenReqs are set to 0 later on.
		minGenerationReq[g] = genPtr->minGenerationReq;
		if (minGenerationReq[g] > min(genPtr->rampUpLim * periodLength, genPtr->rampDownLim * periodLength)) {
			minGenerationReq[g] = min(genPtr->rampUpLim * periodLength, genPtr->rampDownLim * periodLength);
		}
	}

	/* Min Up/Down */
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		minUpTimePeriods[g]		= round(genPtr->minUpTime * 60.0/periodLength);
		minDownTimePeriods[g]	= round(genPtr->minDownTime * 60.0/periodLength);
		
		if (minUpTimePeriods[g] < 1)	minUpTimePeriods[g] = 1;
		if (minDownTimePeriods[g] < 1)	minDownTimePeriods[g] = 1;
	}
	
	/* Select a random permutation of scenarios from the sample space */
	rndPermutation.resize(runParam.numLSScen);
	for (int s=0; s<runParam.numLSScen; s++) {
		rndPermutation[s] = rand() % inst->simulations.vals.size();
//		rndPermutation[s] = s;
	}
	
	
	/* Mean Generator Capacity */
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		auto it = inst->observations["RT"].mapVarNamesToIndex.find(genPtr->name);
		
		bool stocSupply = (it != inst->observations["RT"].mapVarNamesToIndex.end()) ? true : false;
		
		if (stocSupply) {
			/* stochastic supply */
			fill(expCapacity[g].begin(), expCapacity[g].end(), 0);
			
			int period;
			for (int t=0; t<numPeriods; t++) {
				period = (beginMin/periodLength)+(t*numBaseTimePerPeriod);
				
				if ( t==0 && probType==ShortTerm ) {
					it = inst->observations["RT"].mapVarNamesToIndex.find(genPtr->name);
					expCapacity[g][t] = min( inst->observations["RT"].vals[rep][period][it->second], genPtr->maxCapacity );
				}
				else {
					it = inst->simulations.mapVarNamesToIndex.find(genPtr->name);
					for (int s=0; s<runParam.numLSScen; s++) {
						expCapacity[g][t] += min( inst->simulations.vals[ rndPermutation[s] ][period][it->second], genPtr->maxCapacity );
					}
					expCapacity[g][t] *= 1.0/(double)(runParam.numLSScen);
				}
			}
		}
		else {
			/* deterministic supply */
			fill(expCapacity[g].begin(), expCapacity[g].end(), genPtr->maxCapacity);
		}
	}
	 
	
	/* Misc: Generators which are scheduled in ST-UC, must have 0 min-generation-requirement & capacity in DA-UC problem */
	if (probType == DayAhead) {
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst->powSys->generators[g]);
			
			if (!genPtr->isDAUCGen) {
				fill(expCapacity[g].begin(), expCapacity[g].end(), 0.0);
			}
		}
	}
	
	/* Load */
	for (int b=0; b<numBus; b++) {
		Bus *busPtr = &(inst->powSys->buses[b]);
		
		auto it = inst->observations["RT"].mapVarNamesToIndex.find( num2str(busPtr->regionId) );
		
		int period;
		for (int t=0; t<numPeriods; t++) {
			period = (beginMin/periodLength)+(t*numBaseTimePerPeriod);
			
			if ( probType == DayAhead ) {
				busLoad[b][t] = inst->observations["DA"].vals[rep][period][it->second] * busPtr->loadPercentage;
			}
			else {
				busLoad[b][t] = inst->observations["RT"].vals[rep][period][it->second] * busPtr->loadPercentage;
			}
		}
	}
	
	/* Spinning Reserve */
	for (int b=0; b<numBus; b++) {
		for (int t=0; t<numPeriods; t++) {
			busLoad[b][t] *= (1+runParam.spinResPerc);
		}
	}
	
	/* System-level loads */
	fill( sysLoad.begin(), sysLoad.end(), 0.0 );		// reset system load to 0
	for (int t=0; t<numPeriods; t++) {
		for (int b=0; b<numBus; b++) {
			sysLoad[t] += busLoad[b][t];
		}
	}
}


void SUCmaster::formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep) {
	
	/* Initializations */
	this->inst		= &inst;
	this->beginMin	= beginMin;
	this->probType	= probType;
	this->modelType = modelType;
	this->rep		= rep;
	multicut		= false;
	
	/* Prepare Model-Dependent Input Data */
	preprocessing();

	/* Master Formulation */
	// create the variables
	eta = IloNumVarArray (env, multicut ? runParam.numLSScen : 1, 0, IloInfinity, ILOFLOAT);	// 2nd-stage approximation
	
	p = IloArray<IloNumVarArray> (env, numGen);	// production amounts
	IloArray<IloNumVarArray> p_var (env, numGen);	// variable-production amounts

    s = IloArray< IloNumVarArray > (env, numGen);
	x = IloArray< IloNumVarArray > (env, numGen);
	z = IloArray< IloNumVarArray > (env, numGen);
	for (int g=0; g<numGen; g++) {
		s[g] = IloNumVarArray(env, numPeriods, 0, 1, ILOBOOL);
		x[g] = IloNumVarArray(env, numPeriods, 0, 1, ILOBOOL);
		z[g] = IloNumVarArray(env, numPeriods, 0, 1, ILOBOOL);
		p[g]	 = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
		p_var[g] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);

		sprintf(buffer, "s_%d", g);
		s[g].setNames(buffer);
		sprintf(buffer, "x_%d", g);
		x[g].setNames(buffer);
		sprintf(buffer, "z_%d", g);
		z[g].setNames(buffer);
		sprintf(buffer, "p_%d", g);
		p[g].setNames(buffer);
		sprintf(buffer, "pv_%d", g);
		p_var[g].setNames(buffer);

		model.add(s[g]);
		model.add(x[g]);
		model.add(z[g]);
	}
	
	// create the constraints

    // state constraints
    for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		if ( (probType == DayAhead && !genPtr->isDAUCGen) || (probType == ShortTerm && genPtr->isDAUCGen) ) {
			continue;	// skip scheduling constraints for not-to-be-scheduled generators
		}

        // t=0: generators are assumed to be turned on
        model.add( x[g][0] - getGenState(g,-1) == s[g][0] - z[g][0] );
        
        // t>0
        for (int t=1; t<numPeriods; t++) {
            model.add( x[g][t] - x[g][t-1] == s[g][t] - z[g][t] );
        }
    }
    
    // minimum uptime/downtime constraints
    for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		if ( (probType == DayAhead && !genPtr->isDAUCGen) || (probType == ShortTerm && genPtr->isDAUCGen) ) {
			continue;	// skip scheduling constraints for not-to-be-scheduled generators
		}

        // turn on inequalities
        for (int t=1; t<=numPeriods; t++)
        {
            IloExpr lhs (env);
            for (int i = t-minUpTimePeriods[g]+1; i<=t; i++)
            {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += max(0, getGenState(g,i-1) - getGenState(g,i-2));
            }
            model.add( lhs <= x[g][t-1] );
            lhs.end();
        }
    }
    
    for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		if ( (probType == DayAhead && !genPtr->isDAUCGen) || (probType == ShortTerm && genPtr->isDAUCGen) ) {
			continue;	// skip scheduling constraints for not-to-be-scheduled generators
		}

        // turn off inequalities
        for (int t=1; t<=numPeriods; t++)
        {
            IloExpr lhs (env);
            for (int i = t-minDownTimePeriods[g]+1; i<=t; i++)  {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += max(0, getGenState(g,i-1) - getGenState(g,i-2));
            }
            
			if (t-minDownTimePeriods[g]-1 >= 0)	model.add( lhs <= 1 - x[g][t-minDownTimePeriods[g]-1] );
			else								model.add( lhs <= 1 - getGenState(g, (t-minDownTimePeriods[g]-1)));

			lhs.end();
        }
    }
    
    /* demand-based valid inequality *
    // the capacities of operational generators must exceed the system demand at any point
    for (int t=0; t<numPeriods; t++) {
        IloExpr expr (env);
        for (int g=0; g<numGen; g++) {
            expr += expCapacity[g][t] * x[g][t];
        }
        model.add( expr >= sysLoad[t] );
        expr.end();
    }
	/*********************************/
	
	// must-run units must be committed
	for (int g=0; g<numGen; g++) {
		if ( inst.powSys->generators[g].isMustRun ) {
			for (int t=0; t<numPeriods; t++) {
				x[g][t].setBounds(1, 1);
			}
		}
	}
	
	// commit the right set of generators for the right type of problem
	if (probType == DayAhead) {
		/* ST-UC generators will not produce in the DA-UC problem. */
		for (int g=0; g<numGen; g++) {
			if (!(inst.powSys->generators[g].isDAUCGen)) {
				for (int t=0; t<numPeriods; t++) {
					x[g][t].setBounds(0, 0);
				}
			}
		}
	}
	else {
		/* All generators will produce in the ST-UC problem. Commitment
		 * decisions of DA-UC generators will be read from the solution. */
		double genState;
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst.powSys->generators[g]);
			
			if (genPtr->isDAUCGen) {
				for (int t=0; t<numPeriods; t++) {
					genState = getGenState(g, t);
					x[g][t].setBounds(genState, genState);
				}
			}
		}
	}
	
	/****** Symmetry Breaking ******/
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		// compare characteristics of generator g, with an other generator k
		for (int k=g+1; k<numGen; k++) {
			Generator *gen2Ptr = &(inst.powSys->generators[k]);
			
			// check if gen g and k are the same
			if ( genPtr->type == gen2Ptr->type
				&& genPtr->connectedBusName == gen2Ptr->connectedBusName
				&& fabs(genPtr->rampUpLim - gen2Ptr->rampUpLim) < EPSzero
				&& fabs(genPtr->rampDownLim - gen2Ptr->rampDownLim) < EPSzero
				&& fabs(genPtr->minGenerationReq - gen2Ptr->minGenerationReq) < EPSzero
				&& fabs(genPtr->maxCapacity - gen2Ptr->maxCapacity) < EPSzero
				&& fabs(genPtr->minUpTime - gen2Ptr->minUpTime) < EPSzero
				&& fabs(genPtr->minDownTime - gen2Ptr->minDownTime) < EPSzero
				&& fabs(genPtr->isDAUCGen - gen2Ptr->isDAUCGen) < EPSzero
				&& fabs(genPtr->isMustRun - gen2Ptr->isMustRun) < EPSzero
				&& fabs(genPtr->isMustUse - gen2Ptr->isMustUse) < EPSzero
				&& fabs(genPtr->noLoadCost - gen2Ptr->noLoadCost) < EPSzero
				&& fabs(genPtr->startupCost - gen2Ptr->startupCost) < EPSzero
				&& fabs(genPtr->variableCost - gen2Ptr->variableCost) < EPSzero){
				
				// if these generators are the same, create a symmetry-breaking constraint
				for (int t=0; t<numPeriods; t++) {
					model.add( x[g][t] >= x[k][t] );	// if you need to open k, make sure g is already open
				}
				break;
			}
		}
	}
	/****** Symmetry Breaking ******/
	
	// create the objective function
	IloExpr obj (env);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);

		for (int t=0; t<numPeriods; t++) {
			obj += genPtr->startupCost * s[g][t];					// start up cost
			obj += genPtr->noLoadCost*periodLength/60.0 * x[g][t];	// no-load cost
//			obj += genPtr->variableCost*periodLength/60.0 * p[g][t];		// generation cost
			obj += minGenerationReq[g]*genPtr->variableCost*periodLength/60.0 * x[g][t];	// minimum generation cost
		}
	}
	for (int s=0; s<eta.getSize(); s++) {
		obj += 1.0/(double)eta.getSize() * eta[s];
	}
	
	
	/****** MEAN SCENARIO CONSTRAINTS ******/
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			if (genPtr->isMustUse) {
				model.add( p_var[g][t] == (expCapacity[g][t] - minGenerationReq[g]) );
			} else {
				model.add( p_var[g][t] <= (expCapacity[g][t] - minGenerationReq[g]) * x[g][t] );
			}
		}
	}
	
	// production amounts
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++) {
			model.add( p[g][t] == p_var[g][t] + x[g][t] * minGenerationReq[g] );
		}
	}
	
	// ramp up constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		/* t = 0 */
		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {	// prev generation is available
			IloConstraint c ( p[g][t] - getGenProd(g, t-1) <= genPtr->rampUpLim*periodLength );
			sprintf(buffer, "RU_%d_%d", g, t);
			c.setName(buffer);
			model.add(c);
		}
		
		/* t > 0 */
		for (t=1; t<numPeriods; t++) {
			IloConstraint c ( p_var[g][t] - p_var[g][t-1] <= genPtr->rampUpLim*periodLength * x[g][t] - minGenerationReq[g] * s[g][t]);
			sprintf(buffer, "RU_%d_%d", g, t);
			c.setName(buffer);
			model.add(c);
		}
	}	
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		// input-inconsistency check
		int shutDownPeriod = checkShutDownRampDownInconsistency(g);
		
		double rampDownRate;
		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			IloConstraint c ( getGenProd(g, t-1) - p[g][t] <= rampDownRate );
			sprintf(buffer, "RD_%d_%d", g, t); c.setName(buffer); model.add(c);
		}
		for (t=1; t<numPeriods; t++) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			IloConstraint c ( p_var[g][t-1] - p_var[g][t] <= rampDownRate * x[g][t-1] - minGenerationReq[g] * z[g][t] );
			sprintf(buffer, "RD_%d_%d", g, t); c.setName(buffer); model.add(c);
		}
	}

	if (modelType == System)
	{
		IloNumVarArray L (env, numPeriods, 0, IloInfinity, ILOFLOAT);	// load shedding
		IloNumVarArray O (env, numPeriods, 0, IloInfinity, ILOFLOAT);	// over generation
		
		// aggregated-demand constraints
		for (int t=0; t<numPeriods; t++) {
			IloExpr expr (env);
			for (int g=0; g<numGen; g++) {
				expr += p[g][t];
			}
			expr += L[t];
			expr -= O[t];
			model.add( IloRange (env, sysLoad[t], expr, sysLoad[t]) );
		}
		
		// load-shedding penalties
		for (int t=0; t<numPeriods; t++) {
			obj += loadShedPenaltyCoef * L[t];
			obj += overGenPenaltyCoef  * O[t];
		}
	}
	else	// (modelType == Transmission)
	{
		L = IloArray< IloNumVarArray > (env, numBus);	// load shedding
		O = IloArray< IloNumVarArray > (env, numBus);	// over generation
		IloArray< IloNumVarArray > T (env, numBus);		// phase angles
		IloArray< IloNumVarArray > F (env, numLine);	// flows
		
		for (int b=0; b<numBus; b++) {
			L[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
			O[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
			T[b] = IloNumVarArray(env, numPeriods, inst.powSys->buses[b].minPhaseAngle, inst.powSys->buses[b].maxPhaseAngle, ILOFLOAT);
			
			sprintf(buffer, "L_%d", b);
			L[b].setNames(buffer);
			sprintf(buffer, "O_%d", b);
			O[b].setNames(buffer);
			sprintf(buffer, "T_%d", b);
			T[b].setNames(buffer);
		}
		
		for (int l=0; l<numLine; l++) {
			F[l] = IloNumVarArray(env, numPeriods, inst.powSys->lines[l].minFlowLim, inst.powSys->lines[l].maxFlowLim, ILOFLOAT);
			
			sprintf(buffer, "F_%d", l);
			F[l].setNames(buffer);
		}
		
		// No load-shedding in buses with 0 load
		for (int b=0; b<numBus; b++) {
			for (int t=0; t<numPeriods; t++) {
				if (busLoad[b][t] < EPSzero) L[b][t].setUB(0);
			}
		}
		// No over-generation in buses with no generators
		for (int b=0; b<numBus; b++) {
			if (inst.powSys->buses[b].connectedGenerators.size() == 0) {
				for (int t=0; t<numPeriods; t++) {
					O[b][t].setUB(0);
				}
			}
		}
		
		// DC-approximation to AC power flow
		for (int l=0; l<numLine; l++) {
			Line *linePtr = &(inst.powSys->lines[l]);
			
			int orig = linePtr->orig->id;
			int dest = linePtr->dest->id;
			
			for (int t=0; t<numPeriods; t++) {
				model.add( F[l][t] == linePtr->susceptance * (T[orig][t] - T[dest][t]) );
			}
		}
		
		// Flow balance
		for (int b=0; b<numBus; b++) {
			Bus *busPtr = &(inst.powSys->buses[b]);
			
			for (int t=0; t<numPeriods; t++) {
				IloExpr expr (env);
				
				// production (iterate over connected generators)
				for (int g = 0; g < (int) busPtr->connectedGenerators.size(); g++) {
					expr += p[ busPtr->connectedGenerators[g]->id ][t];
				}
				
				// in/out flows (iterate over all arcs)
				for (auto linePtr = inst.powSys->lines.begin(); linePtr != inst.powSys->lines.end(); ++linePtr) {
					if (linePtr->orig->id == busPtr->id) {	// if line-origin and bus have same ids, then this is outgoing flow
						expr -= F[linePtr->id][t];
					}
					if (linePtr->dest->id == busPtr->id) {	// if line-destination and bus have same ids, then this is outgoing flow
						expr += F[linePtr->id][t];
					}
				}
				
				// load shedding
				expr += L[b][t];
				
				// over generation
				expr -= O[b][t];
				
				// constraint
				model.add( expr == busLoad[b][t] );
				
				// free up memory
				expr.end();
			}
		}
		
		IloExpr meanValObjFunction (env);
		for (int t=0; t<numPeriods; t++) {
			for (int g=0; g<numGen; g++) {
				Generator *genPtr = &(inst.powSys->generators[g]);
				meanValObjFunction += genPtr->variableCost*periodLength/60.0 * p_var[g][t];		// variable-generation cost
			}
			for (int b=0; b<numBus; b++) {
				meanValObjFunction += overGenPenaltyCoef * O[b][t];
				meanValObjFunction += loadShedPenaltyCoef * L[b][t];
			}
		}
		model.add(eta[0] >= meanValObjFunction);
	}
	 /**************************/
	
    // set the objective function
	model.add( IloMinimize(env, obj) );
	
    // formulate the subproblem
	recourse.formulate(inst, probType, modelType, beginMin, rep, xvals, rndPermutation, expCapacity);
	
	// formulate the warm-up problem
	warmUpProb.formulate(inst, probType, modelType, beginMin, rep);
	warmUpProb.cplex.setParam(IloCplex::SolnPoolGap, 5e-2);
	warmUpProb.cplex.setParam(IloCplex::SolnPoolCapacity, 5);

    // prepare the solver
	cplex.extract(model);

    /*************************************************************************
     * DISABLED: Didn't improve the progress much.. reduced it in some cases.
     *
     * // assign priorities to state variables
     * for (int g=0; g<numGen; g++) {
	 *    for (int t=0; t<numPeriods; t++) {
     *         cplex.setPriority(x[g][t], 1);
     *     }
     * }
    /************************************************************************/
	
	cplex.setOut( inst.out() );
	cplex.setWarning( inst.out() );
	cplex.setParam(IloCplex::Threads, LShapedMasterCPXThreads);
//	cplex.setParam(IloCplex::ParallelMode, IloCplex::Opportunistic);
//	cplex.setParam(IloCplex::MIPEmphasis, IloCplex::MIPEmphasisBestBound);
    cplex.setParam(IloCplex::EpGap, 1e-2);
}

bool SUCmaster::solve () {
	try{
		bool status;
		
		// Benders' decomposition
		inst->out() << "Executing Benders' decomposition for ";
		if (probType == DayAhead)	inst->out() << "DAUC,";
		else						inst->out() << "STUC,";
		inst->out() << " at time = " << beginMin << " (mins)." << endl;
		
		/****** Process the LP Relaxation ******
		// convert variables from BOOL to FLOAT
		IloNumVarArray vars (env);
		for (int g=0; g<numGen; g++) {
			for (int t=0; t<numPeriods; t++) {
				vars.add(x[g][t]);
				vars.add(s[g][t]);
				vars.add(z[g][t]);
			}
		}
		IloConversion convertToLP (env, vars, ILOFLOAT);
		vars.end();
		
		// add conversion to the model
		model.add(convertToLP);
		
		// add a dummy binary to keep the problem as an "MIP".
		model.add( IloNumVar(env, 0, IloInfinity, ILOBOOL) );
		
		// initialize the LP relaxation
		BendersCuts = IloRangeArray (env);
		LinProgRelaxFlag = true;
		cplex.setParam(IloCplex::TiLim, (probType == DayAhead)*600 + (probType == ShortTerm)*300);
		// solve
		status = cplex.solve();
		
		// reformulate the problem with the new cuts
		model.remove(convertToLP);
		model.add(BendersCuts);
		cplex.extract(model);
		LinProgRelaxFlag = false;
		/****** Process the LP Relaxation ******/
		
		/****** Mean Value Problem ******
		inst->out() << "****** SOLVING MEAN VALUE PROBLEM ******" << endl;
		
		// convert variables from BOOL to FLOAT
		IloNumVarArray vars (env);
		for (int g=0; g<numGen; g++) {
			for (int t=0; t<numPeriods; t++) {
				vars.add(x[g][t]);
				vars.add(s[g][t]);
				vars.add(z[g][t]);
			}
		}
		IloConversion convertToLP (env, vars, ILOFLOAT);
		vars.end();
		
		// add conversion to the model
		model.add(convertToLP);
		
		// add a dummy binary to keep the problem as an "MIP".
		model.add( IloNumVar(env, 0, IloInfinity, ILOBOOL) );
		

		
		MeanProbFlag = true;
		BendersCuts = IloRangeArray (env);
		cplex.solve();
		inst->out() << "Mean-Value Obj = " << cplex.getObjValue() << endl;
		model.remove(convertToLP);
		model.add(BendersCuts);
		cplex.extract(model);
		evaluatedSolns.clear();
		MeanProbFlag = false;			
		/****** Mean Value Problem ******/
		
		LinProgRelaxFlag = false;
		
		/* warm up */
		inst->out() << "Executing warm up MIP..." << endl;
		status = warmUpProb.solve();
		if (status) {
			setWarmUp();
			inst->out() << "Warm up model provided " << warmUpProb.cplex.getSolnPoolNsolns() << " solutions (Best Obj= " << warmUpProb.getObjValue() << ")." << endl;
		} else {
			inst->out() << "Warm up has failed." << endl;
		}
		/*        */

		/****** Process the MIP ******/
		inst->out() << "****** SOLVING THE PROBLEM ******" << endl;
		cplex.setParam(IloCplex::TiLim, (probType == DayAhead)*7200 + (probType == ShortTerm)*1800);
//		cplex.use(LazySepCallback(env, *this));
//		cplex.use(rounding(env, x, *this));
//		if (probType == DayAhead) {
//			cplex.use(IncCallback(env, *this));
//		}
		
		// create an LShapedCallback
		LShapedCallback callback (*this);
		CPXLONG contextMask = 0;
		contextMask |= IloCplex::Callback::Context::Id::Candidate;
		cplex.use(&callback, contextMask);
		
		status = cplex.solve();
		
		if (status) {
			inst->out() << "Optimization is completed with status " << cplex.getCplexStatus() << endl;
			inst->out() << "Obj = \t" << cplex.getObjValue() << endl;
			inst->out() << "LB = \t" << cplex.getBestObjValue() << endl;
		} else {
			cplex.exportModel("infeasible.lp");
			inst->printSolution("infeasible");
			inst->out() << "Benders' decomposition has failed." << endl;
			exit(1);
		}
		/****** Process the MIP ******/
		
		// Record the optimal solution
		if (status) {
			for (int g=0; g<numGen; g++) {
				Generator *genPtr = &(inst->powSys->generators[g]);
				
				for (int t=0; t<numPeriods; t++) {
					if ( (probType == DayAhead && genPtr->isDAUCGen) || (probType == ShortTerm && !genPtr->isDAUCGen) ) {
						setGenState(g,t, cplex.getValue(x[g][t]));
					}
					// Note: "Expected" production amounts are recorded in the callback functions. The below code
					// sets production to 0 when generators are not committed. This is essentially correcting potential
					// numerical errors.
					if ( probType == DayAhead || (probType == ShortTerm && beginMin == 0) ) {
						if (!getGenState(g,t)) setUCGenProd(g, t, 0.0);
					}
				}
			}
		}
		
		return status;
	}
	catch (IloException &e) {
		cout << e << endl;
		return false;
	}
}

void SUCmaster::setWarmUp ()
{
	for (int sol=0; sol<warmUpProb.cplex.getSolnPoolNsolns(); sol++) {
		IloNumVarArray	vars (env);
		IloNumArray		vals (env);

		for (int g=0; g<numGen; g++) {
			for (int t=0; t<numPeriods; t++) {
				vars.add( s[g][t] );
				vals.add( warmUpProb.cplex.getValue(warmUpProb.s[g][t], sol) );
				
				vars.add( x[g][t] );
				vals.add( warmUpProb.cplex.getValue(warmUpProb.x[g][t], sol) );
				
				vars.add( z[g][t] );
				vals.add( warmUpProb.cplex.getValue(warmUpProb.z[g][t], sol) );
			}
		}
		cplex.addMIPStart(vars, vals, IloCplex::MIPStartSolveFixed);
		
		vars.end();
		vals.end();
	}
}

/****************************************************************************
 * getGenState
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded status of the generator.
 ****************************************************************************/
bool SUCmaster::getGenState(int genId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// return the requested generator state
	if (reqSolnComp < 0) {										// all generators are assumed to be online, for a long time, at t=0.
		return true;
	}
	else if (reqSolnComp < (int) inst->solution.x[genId].size()) {	// return the corresponding solution
		return (inst->solution.x[genId][reqSolnComp] > EPSzero);
	}
	else {														// asking what's beyond the planning horizon, we return the last solution
		return (inst->solution.x[genId][ inst->solution.x[genId].size()-1 ] > EPSzero);
	}
}

/****************************************************************************
 * setGenState
 * - Fills the (genId, correspondingComponent) of the Solution.x object.
 ****************************************************************************/
void SUCmaster::setGenState(int genId, int period, double value) {
	// which Solution component is being set?
	int solnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// correct potential numerical errors (important for binary variables)
	value = round(value);

	// set the solution
	if (solnComp >= 0 && solnComp < (int) inst->solution.x[genId].size()) {
		for (int t=solnComp; t<solnComp+numBaseTimePerPeriod; t++) {
			inst->solution.x[genId][t] = value;
		}
	}
	else {
		// Setting generator state at time that is beyond the planning horizon
	}
}

/****************************************************************************
 * getObjValue
 * - Returns the objective value of the last optimization.
 ****************************************************************************/
double SUCmaster::getObjValue() {
	return cplex.getObjValue();
}

/****************************************************************************
 * getEDGenProd
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded generation of the generator by the ED model.
 ****************************************************************************/
double SUCmaster::getEDGenProd(int genId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// return the requested generator state
	if (reqSolnComp < 0) {										// all generators are assumed to be online, for a long time, at t=0.
		cout << "Error: Initial production levels are not available" << endl;
		exit(1);
	}
	else if (reqSolnComp < (int) inst->solution.x[genId].size()) {	// return the corresponding solution
		return inst->solution.g_ED[genId][reqSolnComp];
	}
	else {														// asking what's beyond the planning horizon, we return the last solution
		cout << "Error: Production levels beyond the planning horizon are not available" << endl;
		exit(1);
	}
}

/****************************************************************************
 * setUCGenProd
 * - Fills the (genId, correspondingComponent) of the Solution.g_UC object.
 ****************************************************************************/
void SUCmaster::setUCGenProd(int genId, int period, double value) {
	// which Solution component is being set?
	int solnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// correct potential numerical errors
	value = max(0.0, value);

	// set the solution
	if (solnComp >= 0 && solnComp < (int) inst->solution.g_UC[genId].size()) {
		for (int t=solnComp; t<solnComp+numBaseTimePerPeriod; t++) {
			inst->solution.g_UC[genId][t] = value;
		}
	}
	else {
		// Setting generator production at time that is beyond the planning horizon
	}
}

/****************************************************************************
 * getGenProd
 * - If period != the beginning of the planning horizon, returns getEDGenProd
 * - If period is the beginning of the planning horizon:
 * 		- Use generator histories from the previous run:
 *		- Don't use generator histories
 *			- DA-UC: returns -infinity due to lack of info
 *			- ST-UC: returns 1st-period generation amount from the DA solve,
 * 			if the gen is DA, -infinity if it is ST.
 ****************************************************************************/
double SUCmaster::getGenProd(int g, int t) {
	if (t < 0 && beginMin == 0) {
		if ( runParam.useGenHistory ) {
			cout << "Not implemented yet" << endl;
		}
		else {
			Generator *genPtr = &(inst->powSys->generators[g]);
			if (probType == ShortTerm && genPtr->isDAUCGen) {
				return getUCGenProd(g, 0);
			}
		}
	} else {
		return getEDGenProd(g, t);
	}
	return -INFINITY;
}

/****************************************************************************
 * checkShutDownRampDownInconsistency
 * This function assures the consistency of the model inputs with the model's
 * feasible set. It checks if generator g is not be able to shut down within
 * the time horizon of the model, due to a high starting-generation-amount
 * and tight ramping-down restrictions.
 * - If an inconsistency is detected, the shut down period is returned so
 * that ramping restrictions are relaxed at the time of the shut down.
 * - The function returns -1 if it cannot detect an inconsistency.
 
 ****************************************************************************/
int SUCmaster::checkShutDownRampDownInconsistency (int g) {
	// This is a problem pertaining only to DA generators in the ST-UC problem
	if (probType==DayAhead || !inst->powSys->generators[g].isDAUCGen) return -1;
	
	// Was g producing at period -1?
	if ( getGenProd(g, -1) <= EPSzero )	return -1;
	
	// Will g shutdown during the problem's planning horizon?
	bool answer = false;
	int t;
	for (t=0; t<numPeriods; t++) {
		if (!getGenState(g, t)) {
			answer = true;
			break;
		}
	}
	if (!answer) return -1;
	
	// Can the generator ramp down?
	if ( inst->powSys->generators[g].rampDownLim * periodLength * (double)(t+1) < getGenProd(g, -1) ) {
		inst->out() << "Warning: Generator " << g << " (" << inst->powSys->generators[g].name << ") cannot ramp down to 0 in the ST-UC problem" << endl;
		return t;
	}
	else {
		return -1;
	}
}

/****************************************************************************
 * getUCGenProd
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded generation of the generator by the DA model.
 ****************************************************************************/
double SUCmaster::getUCGenProd(int genId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// return the requested generator state
	if (reqSolnComp < 0) {										// all generators are assumed to be online, for a long time, at t=0.
		cout << "Error: Initial production levels are not available" << endl;
		exit(1);
	}
	else if (reqSolnComp < (int) inst->solution.x[genId].size()) {	// return the corresponding solution
		return inst->solution.g_UC[genId][reqSolnComp];
	}
	else {														// asking what's beyond the planning horizon, we return the last solution
		cout << "Error: Production levels beyond the planning horizon are not available" << endl;
		exit(1);
	}
}
