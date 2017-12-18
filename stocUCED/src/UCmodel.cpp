//
//  UCmodel.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 11/11/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "UCmodel.hpp"

//TODO: Why is this runParam here?
extern runType runParam;

UCmodel::UCmodel () {
	model = IloModel(env);
	cplex = IloCplex(env);
}

UCmodel::~UCmodel() {
	env.end();
}

/****************************************************************************
 * preprocessing
 * - Initializes basic parameters
 * - Fills the following containers, according to model units & assumptions:
 *	 minGenerationReq
 *   minUpTimePeriods
 *   minDownTimePeriods
 *   expCapacity
 * - If it is a DA-UC problem, sets the capacity and minGenerationReq to 0
 * for generators which must be scheduled at ST-UC.
 *
 * - Assumption 1: Min generation requirement must be less than ramping
 * rates, otherwise, generators cannot switch on/off.
 * - Assumption 2: Min up/down times must at least be 1 period.
 ****************************************************************************/
void UCmodel::preprocessing ()
{
	// initialize basic parameters
	numGen	   = inst->powSys->numGen;
	numLine    = inst->powSys->numLine;
	numBus     = inst->powSys->numBus;
	
	numPeriods	 = (probType == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
	periodLength = (probType == DayAhead) ? runParam.DA_resolution : runParam.ST_resolution;
	
	numSolnCompsPerPeriod = (int)round(periodLength) / runParam.ED_resolution;

	
	// initialize the containers
	minGenerationReq.resize(inst->powSys->numGen);		// minimum production requirements
	minUpTimePeriods.resize(inst->powSys->numGen);		// minimum uptime in periods
	minDownTimePeriods.resize(inst->powSys->numGen);	// minimum downtime in periods
	resize_matrix(expCapacity, numGen, numPeriods);		// mean generator capacities

	// minGenerationReq & Assumption 1
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		minGenerationReq[g] = genPtr->minGenerationReq * periodLength/60.0;
		if (minGenerationReq[g] > min(genPtr->rampUpLim * periodLength, genPtr->rampDownLim * periodLength)) {
			minGenerationReq[g] = min(genPtr->rampUpLim * periodLength, genPtr->rampDownLim * periodLength);
		}
	}
	
	// minUpTimePeriods, minDownTimePeriods & Assumption 2
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		minUpTimePeriods[g]		= round(genPtr->minUpTime * 60.0/periodLength);
		minDownTimePeriods[g]	= round(genPtr->minDownTime * 60.0/periodLength);

		if (minUpTimePeriods[g] < 1)	minUpTimePeriods[g] = 1;
		if (minDownTimePeriods[g] < 1)	minDownTimePeriods[g] = 1;
	}
	
	// expected capacity
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		if (genPtr->type == Generator::WIND || genPtr->type == Generator::SOLAR) {
			/* random supply */
			// TODO: Get the expectation of the random supply
			fill(expCapacity[g].begin(), expCapacity[g].end(), genPtr->maxCapacity*periodLength/60.0);
		}
		else {
			/* deterministic supply */
			fill(expCapacity[g].begin(), expCapacity[g].end(), genPtr->maxCapacity*periodLength/60.0);
		}
	}

	// Generators which are scheduled in ST-UC, must have 0 min-generation-requirement & capacity in DA-UC problem
	if (probType == DayAhead) {
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst->powSys->generators[g]);
			
			if (!genPtr->isBaseLoadGen) {
				minGenerationReq[g] = 0.0;
				fill(expCapacity[g].begin(), expCapacity[g].end(), 0.0);
			}
		}
	}
}

	//	aggregated_demand.resize(nb_periods);
//
//	resize_matrix(demand,   inst->numBus, nb_periods);
//	resize_matrix(capacity, inst->numGen, nb_periods);
//
//	// set basic parameters
//	int begin_period = begin_hour * 60 / period_len;
//
//	/* Calculate demand for each 'period' and at each 'bus'
//	 *
//	for (int b=0; b<inst->numBus; b++) {
//		for (int t=0; t<nb_periods; t++)
//		{
//			if ( prob_type == DayAhead) {
//				demand[b][t] = inst->load_perHour   [ inst->bus_regionId[b] ][begin_period+t] * inst->bus_loadPerc[b];
//			} else {
//				demand[b][t] = inst->load_perSubHour[ inst->bus_regionId[b] ][begin_period+t] * inst->bus_loadPerc[b];
//			}
//		}
//	}
//
//	// Calculate aggregated demand for each period
//	fill( aggregated_demand.begin(), aggregated_demand.end(), 0.0 );	// initialize to 0
//	for (int r=0; r<inst->load_perHour.size(); r++) {
//		for (int t=0; t<nb_periods; t++) {
//			if ( prob_type == DayAhead ) {
//				aggregated_demand[t] += inst->load_perHour   [r][begin_period+t];
//			} else {
//				aggregated_demand[t] += inst->load_perSubHour[r][begin_period+t];
//			}
//		}
//	}

/****************************************************************************
 * formulate
 * - Performs preprocessing on the instance, to compute necessary model
 * parameters.
 * - Formulates the mathematical program.
 ****************************************************************************/
void UCmodel::formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin) {
	
	/* Initializations */
	this->inst		= &inst;
	this->beginMin	= beginMin;
	this->probType	= probType;


	/* Prepare Model-Dependent Input Data */
	preprocessing();
	
	
	/* Prepare the Mathematical Model */
	
	/* Decision Variables */
	s	  = IloArray<IloNumVarArray> (env, numGen);	// start up vars
	x	  = IloArray<IloNumVarArray> (env, numGen);	// state vars
	z	  = IloArray<IloNumVarArray> (env, numGen);	// shut down vars
	p	  = IloArray<IloNumVarArray> (env, numGen);	// production amounts
	p_var = IloArray<IloNumVarArray> (env, numGen);	// variable-production amounts
	
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


	/**** Constraints (Traditional Formulation) ****/
	
	// state constraints
	for (int g=0; g<numGen; g++) {
		// t=0: generators are assumed to be turned on
		model.add( x[g][0] - getGenState(g,-1) == s[g][0] - z[g][0] );
	
		// t>0
		for (int t=1; t<numPeriods; t++) {
			model.add( x[g][t] - x[g][t-1] == s[g][t] - z[g][t] );
		}
	}
	
	// minimum uptime constraints
	for (int g=0; g<numGen; g++) {
		// turn on inequalities
		for (int t=1; t<=numPeriods; t++)
		{
			IloExpr lhs (env);
			for (int i = t-minUpTimePeriods[g]+1; i<=t; i++) {
				if (i-1 >= 0)	lhs += s[g][i-1];
				else			lhs += max(0, getGenState(g,i-1) - getGenState(g,i-2));
			}
			model.add( lhs <= x[g][t-1] );
			lhs.end();
		}
	}
	
	// minimum downtime constraints
	for (int g=0; g<numGen; g++) {
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
	
	// must-run units must be committed
	for (int g=0; g<numGen; g++) {
		if ( inst.powSys->generators[g].isMustRun ) {
			for (int t=0; t<numPeriods; t++) {
				x[g][t].setLB(1);
			}
		}
	}
	
	// commit the right set of generators for the right type of problem
	if (probType == DayAhead) {
		/* ST-UC generators will not produce in the DA-UC problem. */
		for (int g=0; g<numGen; g++) {
			if (!(inst.powSys->generators[g].isBaseLoadGen)) {
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
			
			if (genPtr->isBaseLoadGen) {
				for (int t=0; t<numPeriods; t++) {
					genState = getGenState(g, t);
					x[g][t].setBounds(genState, genState);
				}
			}
		}
	}
	
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			if (genPtr->isMustRun) {
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
		
		for (int t=1; t<numPeriods; t++) {
			model.add( p_var[g][t] - p_var[g][t-1] <= genPtr->rampUpLim*periodLength * x[g][t] - minGenerationReq[g] * s[g][t]);
		}
	}
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++) {
			model.add( p_var[g][t-1] - p_var[g][t] <= genPtr->rampDownLim*periodLength * x[g][t-1] - minGenerationReq[g] * z[g][t] );
		}
	}
	
	
	/** Objective Function **/
	IloExpr obj (env);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			obj += genPtr->startupCost * s[g][t];							// start up cost
			obj += minGenerationReq[g] * genPtr->variableCost * x[g][t];	// cost of producing minimum production amount
			obj += genPtr->noLoadCost*periodLength/60.0 * x[g][t];			// no-load cost
			obj += genPtr->variableCost * p_var[g][t];						// variable cost
		}
	}
	
	
	/** Model-dependent obj function components, constraints **/
	if (modelType == System)
	{
		IloNumVarArray L (env, numPeriods, 0, IloInfinity, ILOFLOAT);	// load shedding
		
		// aggregated-demand constraints
		for (int t=0; t<numPeriods; t++) {
			IloExpr expr (env);
			for (int g=0; g<numGen; g++) {
				expr += p[g][t];
			}
			expr += L[t];
			model.add( IloRange (env, aggregated_demand[t], expr) );
		}
		
		// load-shedding penalties
		for (int t=0; t<numPeriods; t++) {
			obj += loadShedPenaltyCoef * L[t];
		}
	}
	else	// (modelType == Transmission)
	{
		IloArray< IloNumVarArray > L (env, numBus);		// load shedding
		IloArray< IloNumVarArray > T (env, numBus);		// phase angles
		IloArray< IloNumVarArray > F (env, numLine);	// flows

		for (int b=0; b<numBus; b++) {
			L[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
			T[b] = IloNumVarArray(env, numPeriods, 0, 2*pi, ILOFLOAT);
			
			sprintf(buffer, "L_%d", b);
			L[b].setNames(buffer);
			sprintf(buffer, "T_%d", b);
			T[b].setNames(buffer);
		}
		
		for (int l=0; l<numLine; l++) {
			F[l] = IloNumVarArray(env, numPeriods,
								  inst.powSys->lines[l].minFlowLim*periodLength/60.0,
								  inst.powSys->lines[l].maxFlowLim*periodLength/60.0, ILOFLOAT);
			
			sprintf(buffer, "F_%d", l);
			F[l].setNames(buffer);
		}
		
		// TODO: No Load-shedding in buses with 0 load
		
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
				for (int g=0; g<busPtr->connectedGenerators.size(); g++) {
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
				
				// constraint
				//TODO: set the demands
				model.add( expr == 0.0 );
				
				// free up memory
				expr.end();
			}
		}
		
		// load-shedding penalties
		for (int b=0; b<numBus; b++) {
			for (int t=0; t<numPeriods; t++) {
				obj += loadShedPenaltyCoef * L[b][t];
			}
		}
	}
	
	/** Finalize **/
	model.add( IloMinimize(env, obj) );
	obj.end();
	
	cplex.extract(model);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::EpGap, 1e-2);
}

/****************************************************************************
 * solve
 * - Solves the model and records the solution into inst->solution.
 ****************************************************************************/
bool UCmodel::solve() {
	
	bool status;
	try {
		status = cplex.solve();
		
		// record the solution
		if (status) {
			for (int g=0; g<numGen; g++) {
				Generator *genPtr = &(inst->powSys->generators[g]);
				
				if ( (probType == DayAhead && genPtr->isBaseLoadGen) || (probType == ShortTerm && !genPtr->isBaseLoadGen) ) {
					for (int t=0; t<numPeriods; t++) {
						setGenState(g,t, cplex.getValue(x[g][t]));
					}
				}
			}
		}
	}
	catch (IloException &e) {
		cout << e << endl;
	}
	
	return status;
}

//bool UCmodel::solve () {
//	try{
//		bool status = cplex.solve();
//
//		cout << "Optimization is completed with status " << cplex.getCplexStatus() << endl;
//
//		// if the problem has a solution, record it
//		if (status) {
//			soln.allocateMem(inst->numGen, nb_periods);
//
//			for (int g=0; g<inst->numGen; g++) {
//				for (int t=0; t<nb_periods; t++) {
//					soln.s[g][t] = cplex.getValue(s[g][t]);
//					soln.x[g][t] = cplex.getValue(x[g][t]);
//					soln.z[g][t] = cplex.getValue(z[g][t]);
//				}
//			}
//
//			solnPool.resize( cplex.getSolnPoolNsolns() );
//			for (int sol=0; sol<cplex.getSolnPoolNsolns(); sol++) {
//				solnPool[sol].allocateMem(inst->numGen, nb_periods);
//
//				for (int g=0; g<inst->numGen; g++) {
//					for (int t=0; t<nb_periods; t++) {
//						solnPool[sol].s[g][t] = cplex.getValue(s[g][t], sol);
//						solnPool[sol].x[g][t] = cplex.getValue(x[g][t], sol);
//						solnPool[sol].z[g][t] = cplex.getValue(z[g][t], sol);
//					}
//				}
//			}
//
//			// update generator status for rolling horizon
//			inst->updateGenStatus(prob_type, soln.x);
//		}
//
//		return status;
//	}
//	catch (IloException &e) {
//		cout << e << endl;
//		return false;
//	}
//}

/****************************************************************************
 * getGenState
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded status of the generator.
 ****************************************************************************/
bool UCmodel::getGenState(int genId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numSolnCompsPerPeriod;
	
	// return the requested generator state
	if (reqSolnComp < 0) {										// all generators are assumed to be online, for a long time, at t=0.
		return true;
	}
	else if (reqSolnComp < inst->solution.x[genId].size()) {	// return the corresponding solution
		return round(inst->solution.x[genId][reqSolnComp]);
	}
	else {														// error
		cout << "Error: You cannot access a solution component that is beyond the planning horizon" << endl;
		exit(-1);
	}
}

/****************************************************************************
 * setGenState
 * - Fills the (genId, correspondingComponent) of the Solution.x object.
 ****************************************************************************/
void UCmodel::setGenState(int genId, int period, double value) {
	// which Solution component is being set?
	int solnComp = beginMin/runParam.ED_resolution + period*numSolnCompsPerPeriod;
	
	// set the solution
	if (solnComp >= 0 && solnComp < inst->solution.x[genId].size()) {
		for (int t=solnComp; t<solnComp+numSolnCompsPerPeriod; t++) {
			inst->solution.x[genId][t] = value;
		}
	}
	else {
		cout << "Error: Setting generator state out of the bounds of the horizon" << endl;
		exit(-1);
	}
}
