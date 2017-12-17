//
//  UCmodel.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 11/11/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "UCmodel.hpp"

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
 * - Fills the following containers, according to model units & assumptions:
 *	 minGenerationReq
 *   minUpTimePeriods
 *   minDownTimePeriods
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
	numPeriods = 24;	//TODO: clean this up
	
	// initialize the containers
	minGenerationReq.resize(inst->powSys->numGen);		// minimum production requirements
	minUpTimePeriods.resize(inst->powSys->numGen);		// minimum uptime in periods
	minDownTimePeriods.resize(inst->powSys->numGen);	// minimum downtime in periods

	// minGenerationReq & Assumption 1
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		minGenerationReq[g] = genPtr->minGenerationReq * (double)period_len/60.0;
		if (minGenerationReq[g] > min(genPtr->rampUpLim * (double)period_len, genPtr->rampDownLim * (double)period_len)) {
			minGenerationReq[g] = min(genPtr->rampUpLim * (double)period_len, genPtr->rampDownLim * (double)period_len);
		}
	}
	
	// minUpTimePeriods, minDownTimePeriods & Assumption 2
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		minUpTimePeriods[g]		= round(genPtr->minUpTime * 60.0/(double)period_len);
		minDownTimePeriods[g]	= round(genPtr->minDownTime * 60.0/(double)period_len);

		if (minUpTimePeriods[g] < 1)	minUpTimePeriods[g] = 1;
		if (minDownTimePeriods[g] < 1)	minDownTimePeriods[g] = 1;
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
//
//	// Calculate hourly generator capacities, take the average of random-supply realizations
//	for (int g = 0; g<inst->numGen; g++) {
//		if ( inst->rndSupply_gen[g] ) {	// take the average of the supply realizations
//			for (int t = 0; t<nb_periods; t++) {
//				capacity[g][t] = inst->sceProb[0] * inst->rndSupply[g][0][begin_period+t] * (double)period_len/60.0;
//				for (int s = 1; s < inst->numScen; s++) {
//					capacity[g][t] += inst->sceProb[s] * inst->rndSupply[g][s][begin_period+t] * (double)period_len/60.0;
//				}
//			}
//		}
//		else {							// use the capacities
//			fill ( capacity[g].begin(), capacity[g].end(), inst->capacity[g] * (double)period_len/60.0 );
//		}
//	}
//
//	// If this is a DA-UC problem, set the ST-UC generators' capacities to 0
//	if (prob_type == DayAhead) {
//		for (int g=0; g<inst->numGen; g++) {
//			if (!inst->dayahead_gen[g]) {
//				minimum_production_req[g] = 0.0;
//				fill( capacity[g].begin(), capacity[g].end(), 0 );
//			}
//		}
//	}
//}

/****************************************************************************
 * formulate
 * - Performs preprocessing on the instance, to compute necessary model
 * parameters.
 * - Formulates the mathematical program.
 ****************************************************************************/
void UCmodel::formulate (instance &inst, ProblemType prob_type, ModelType modelType, int begin_hour) {
	
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
	//	for (int g=0; g<inst.numGen; g++) {
	//
	//		// t=0: generators are assumed to be turned on
	//		model.add( x[g][0] - inst.getGenStatus(g, begin_period-1) == s[g][0] - z[g][0] );
	//
	//		// t>0
	//		for (int t=1; t<nb_periods; t++) {
	//			model.add( x[g][t] - x[g][t-1] == s[g][t] - z[g][t] );
	//		}
	//	}
	//
	//	// minimum uptime constraints
	//	for (int g=0; g<inst.numGen; g++)
	//	{
	//		// turn on inequalities
	//		for (int t=1; t<=nb_periods; t++)
	//		{
	//			IloExpr lhs (env);
	//			for (int i = t-minimum_uptime_periods[g]+1; i<=t; i++)
	//			{
	//				if (i-1 >= 0)	lhs += s[g][i-1];
	//				else			lhs += max(0, inst.getGenStatus(g, begin_period+i-1) - inst.getGenStatus(g, begin_period+i-2));
	//			}
	//			model.add( lhs <= x[g][t-1] );
	//			lhs.end();
	//		}
	//	}
	//
	//	// minimum downtime constraints
	//	for (int g=0; g<inst.numGen; g++)
	//	{
	//		// turn off inequalities
	//		for (int t=1; t<=nb_periods; t++)
	//		{
	//			IloExpr lhs (env);
	//			for (int i = t-minimum_dotime_periods[g]+1; i<=t; i++)  {
	//				if (i-1 >= 0)	lhs += s[g][i-1];
	//				else			lhs += 0;	// otherwise, generator is assumed to be operational (but turned on way earlier in the past)
	//			}
	//
	//			if (t-minimum_dotime_periods[g]-1 >= 0)	model.add( lhs <= 1 - x[g][t-minimum_dotime_periods[g]-1] );
	//			else									model.add( lhs <= 1 - inst.getGenStatus(g, (begin_period+t-minimum_dotime_periods[g]-1)));
	//			lhs.end();
	//		}
	//	}
	
	// must-run units must be committed
	for (int g=0; g<numGen; g++) {
		if ( inst.powSys->generators[g].isMustRun ) {
			for (int t=0; t<numPeriods; t++) {
				x[g][t].setLB(1);
			}
		}
	}
	
	// commit the right set of generators for the right type of problem
	if (prob_type == DayAhead) {
		/* ST-UC generators will not produce in the DA-UC problem. */
		for (int g=0; g<numGen; g++) {
			if (!(inst.powSys->generators[g].isBaseLoadGen)) {
				for (int t=0; t<numPeriods; t++) {
					x[g][t].setBounds(0, 0);
				}
			}
		}
	}
	else
	{
		//TODO:
		/* All generators will produce in the ST-UC problem. Commitment
		 * decisions of DA-UC generators will be read from file.
		 */
	}
	//		ifstream input;
	//		input.open("./DayAhead.sol");
	//		if (input.fail()) {
	//			cout << "Cannot find the DA-UC solution" << endl;
	//		}
	//
	//		vector< vector<double> > DAUC_commitments;
	//		resize_matrix(DAUC_commitments, inst.numGen, runParam.DA_numPeriods);
	//
	//		for (int g=0; g<inst.numGen; g++) {
	//			input >> DAUC_commitments[g][0];	// will be overwritten
	//			for (int t = 0; t < runParam.DA_numPeriods; t++) {
	//				input >> DAUC_commitments[g][t];
	//			}
	//		}
	//		input.close();
	//
	//		// read DA-UC generators' commitment decisions from file
	//		for (int g=0; g<inst.numGen; g++) {
	//			if (inst.dayahead_gen[g]) {
	//				for (int t=0, hour=begin_hour; t<nb_periods; t++, hour=begin_hour+t*period_len/60) {
	//					x[g][t].setBounds( DAUC_commitments[g][hour], DAUC_commitments[g][hour] );
	//				}
	//			}
	//		}
	//	}
	//
	
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			if (genPtr->isMustRun) {
				model.add( p_var[g][t] == (capacity[g][t] - minGenerationReq[g]) );
			} else {
				model.add( p_var[g][t] <= (capacity[g][t] - minGenerationReq[g]) * x[g][t] );
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
			model.add( p_var[g][t] - p_var[g][t-1] <= genPtr->rampUpLim*(double)period_len * x[g][t] - minGenerationReq[g] * s[g][t]);
		}
	}
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++) {
			model.add( p_var[g][t-1] - p_var[g][t] <= genPtr->rampDownLim*(double)period_len * x[g][t-1] - minGenerationReq[g] * z[g][t] );
		}
	}
	
	
	/** Objective Function **/
	IloExpr obj (env);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			obj += genPtr->startupCost * s[g][t];							// start up cost
			obj += minGenerationReq[g] * genPtr->variableCost * x[g][t];	// cost of producing minimum production amount
			obj += genPtr->noLoadCost*(double)period_len/60.0 * x[g][t];	// no-load cost
			obj += genPtr->variableCost * p_var[g][t];						// variable cost
		}
	}
	
	
	/** Model-dependent obj function components, constraints **/
	if (modelType == System) {
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
	else if (modelType == Transmission) {
		IloArray< IloNumVarArray > L (env, numBus);		// load shedding
		IloArray< IloNumVarArray > T (env, numBus);		// phase angles
		IloArray< IloNumVarArray > F (env, numLine);	// flows

		for (int b=0; b<numPeriods; b++) {
			L[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
			T[b] = IloNumVarArray(env, numPeriods, 0, 2*pi, ILOFLOAT);
			
			sprintf(buffer, "L_%d", b);
			L[b].setNames(buffer);
			sprintf(buffer, "T_%d", b);
			T[b].setNames(buffer);
		}
		
		for (int l=0; l<numLine; l++) {
			F[l] = IloNumVarArray(env, numPeriods,
								  inst.powSys->lines[l].minFlowLim*(double)period_len/60.0,
								  inst.powSys->lines[l].maxFlowLim*(double)period_len/60.0, ILOFLOAT);
			
			sprintf(buffer, "F_%d", l);
			F[l].setNames(buffer);
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
		for (auto busPtr = inst.powSys->buses.begin(); busPtr != inst.powSys->buses.end(); ++busPtr) {	// iterate over buses
			for (int t=0; t<numPeriods; t++) {
				
				IloExpr expr (env);
				
				// production (iterate over connected generators)
				for (int g=0; g<busPtr->connectedGenerators.size(); g++) {
					expr += p[ busPtr->connectedGenerators[g]->id ][t];
					cout << "Bus " << busPtr->id << " has " << busPtr->connectedGenerators[g]->name << " connected to it" << endl;
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
				expr += L[busPtr->id][t];
				
				// constraint
				model.add( expr == demand[busPtr->id][t] );
				
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
	else {
		cout << "Invalid model number" << endl;
		exit(-1);
	}
	
	
	/** Finalize **/
	model.add( IloMinimize(env, obj) );
	
	cplex.extract(model);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::EpGap, 1e-2);

}
//	// initialize the instance
//	this->inst	= &inst;
//
//	// initialize model parameters
//	this->begin_hour = begin_hour;
//	this->nb_periods = (prob_type == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
//	this->period_len = (prob_type == DayAhead) ? runParam.DA_resolution : runParam.ST_resolution;
//	this->prob_type  = prob_type;
//
//	int begin_period = begin_hour * 60/period_len;

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
//
//void UCmodel::exportModel()
//{
//	cplex.exportModel("mip.lp");
//}
//
//void UCmodel::printSolution ()
//{
//	char filename[256];
//
//	// TODO: Undo after the instance class is updated.
////	if (prob_type == DayAhead) {
////		sprintf(filename, "%sDayAhead.sol", inst->getName().c_str());
////	} else {
////		sprintf(filename, "%sShortTerm%d.sol", inst->getName().c_str(), begin_hour);
////	}
//
//	ofstream output;
//	output.open(filename);
//
//	for (int g=0; g<inst->numGen; g++) {
//		output << g << "\t";
//		for (int t=0; t<nb_periods; t++) {
//			output << round(soln.x[g][t]) << "\t";
//		}
//		output << endl;
//	}
//
//	output.close ();
//}
//
//void UCmodel::updateSoln (solution &soln)
//{
//	int nb_subhours = 60 / runParam.ST_resolution;
//
//	for (int g=0; g<inst->numGen; g++)
//	{
//		if (prob_type == DayAhead && inst->dayahead_gen[g])	// if the problem was DA-UC, and this generator is a base-load gen, ...
//		{
//			for (int h = 0; h < runParam.DA_numPeriods; h++) {
//				for (int t=0; t<nb_subhours; t++) {
//					soln.x[g][ begin_hour*60/runParam.DA_resolution + h*nb_subhours + t ] = round(this->soln.x[g][h]);
//					soln.s[g][ begin_hour*60/runParam.DA_resolution + h*nb_subhours + t ] = round(this->soln.s[g][h]);
//					soln.z[g][ begin_hour*60/runParam.DA_resolution + h*nb_subhours + t ] = round(this->soln.z[g][h]);
//				}
//			}
//		}
//		else if (prob_type == ShortTerm && !inst->dayahead_gen[g])
//		{
//			for (int t = 0; t < runParam.ST_numPeriods; t++) {
//				soln.x[g][ begin_hour*nb_subhours + t ] = round(this->soln.x[g][t]);
//				soln.s[g][ begin_hour*nb_subhours + t ] = round(this->soln.s[g][t]);
//				soln.z[g][ begin_hour*nb_subhours + t ] = round(this->soln.z[g][t]);
//			}
//		}
//	}
//}
//*/
