////
////  UCmodel.cpp
////  Stoch-UC
////
////  Created by Semih Atakan on 11/11/17.
////  Copyright © 2017 University of Southern California. All rights reserved.
////
///*
//#include "UCmodel.hpp"
//
//extern runType runParam;
//
//UCmodel::UCmodel () {
//	model = IloModel(env);
//	cplex = IloCplex(env);
//}
//
//UCmodel::~UCmodel() {
//	env.end();
//}
//
//void UCmodel::preprocessing ()
//{
//	// initialize the containers
//	minimum_production_req.resize(inst->numGen);
//	minimum_uptime_periods.resize(inst->numGen);
//	minimum_dotime_periods.resize(inst->numGen);
//
//	aggregated_demand.resize(nb_periods);
//
//	resize_matrix(demand,   inst->numBus, nb_periods);
//	resize_matrix(capacity, inst->numGen, nb_periods);
//
//	// set basic parameters
//	int begin_period = begin_hour * 60 / period_len;
//
//	/* Assumption 1: min prod limit must be at most the ramping rates,
//	 * otherwise, generators cannot switch on/off status.
//	 *
//	for (int g=0; g<inst->numGen; g++) {
//		minimum_production_req[g] = inst->min_prod_lim[g] * (double)period_len / 60.0;
//
//		if (minimum_production_req[g] > min(inst->ramp_u_lim[g] * (double)period_len, inst->ramp_d_lim[g] * (double)period_len))	{
//			minimum_production_req[g] = min(inst->ramp_u_lim[g] * (double)period_len, inst->ramp_d_lim[g] * (double)period_len);
//		}
//	}
//
//	/* Assumption 2-3: Minimum up or down times must be at least 1 period
//	 * Note that this may not be guaranteed by 'healthy' input data when we work with subhours
//	 *
//	for (int g=0; g<inst->numGen; g++) {
//		minimum_uptime_periods[g] = round(inst->min_u_time[g] * 60.0 / (double)period_len);
//		minimum_dotime_periods[g] = round(inst->min_d_time[g] * 60.0 / (double)period_len);
//
//		if (minimum_uptime_periods[g] < 1) minimum_uptime_periods[g] = 1;
//		if (minimum_dotime_periods[g] < 1) minimum_dotime_periods[g] = 1;
//	}
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
//
//void UCmodel::formulate (instance &inst, ProblemType prob_type, ModelType model_type, int begin_hour) {
//
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
//
//	/* Prepare Input Data *
//	preprocessing();
//
//	/* Prepare the Mathematical Model *
//
//	/* Decision Variables *
//	s	  = IloArray< IloNumVarArray > (env, inst.numGen);	// start up vars
//	x	  = IloArray< IloNumVarArray > (env, inst.numGen);	// state vars
//	z	  = IloArray< IloNumVarArray > (env, inst.numGen);	// shut down vars
//	p	  = IloArray< IloNumVarArray > (env, inst.numGen);	// production amounts
//	p_var = IloArray< IloNumVarArray > (env, inst.numGen);	// variable-production amounts
//
//	for (int g=0; g<inst.numGen; g++) {
//		s[g] = IloNumVarArray(env, nb_periods, 0, 1, ILOBOOL);
//		x[g] = IloNumVarArray(env, nb_periods, 0, 1, ILOBOOL);
//		z[g] = IloNumVarArray(env, nb_periods, 0, 1, ILOBOOL);
//
//		p[g]	 = IloNumVarArray(env, nb_periods, 0, IloInfinity, ILOFLOAT);
//		p_var[g] = IloNumVarArray(env, nb_periods, 0, IloInfinity, ILOFLOAT);
//
//		char buffer[30];
//
//		sprintf(buffer, "s_%d", g);
//		s[g].setNames(buffer);
//		sprintf(buffer, "x_%d", g);
//		x[g].setNames(buffer);
//		sprintf(buffer, "z_%d", g);
//		z[g].setNames(buffer);
//		sprintf(buffer, "p_%d", g);
//		p[g].setNames(buffer);
//		sprintf(buffer, "pv_%d", g);
//		p_var[g].setNames(buffer);
//
//		model.add(s[g]);
//		model.add(x[g]);
//		model.add(z[g]);
//	}
//
//
//	/**** Constraints (Traditional Formulation) ****
//	// state constraints
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
//
//	// must-run units must be committed
//	for (int g=0; g<inst.numGen; g++) {
//		if (inst.must_run[g]) {
//			for (int t=0; t<nb_periods; t++) {
//				x[g][t].setLB(1);
//			}
//		}
//	}
//
//	// commit the right set of generators for the right type of problem
//	if (prob_type == DayAhead) {
//		/* ST-UC generators will not produce in the DA-UC problem. *
//		for (int g=0; g<inst.numGen; g++) {
//			if (!inst.dayahead_gen[g]) {
//				for (int t=0; t<nb_periods; t++) {
//					x[g][t].setBounds(0, 0);
//				}
//			}
//		}
//	}
//	else
//	{
//		/* All generators will produce in the ST-UC problem. Commitment
//		 * decisions of DA-UC generators will be read from file.
//		 *
//
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
//
//	// capacity constraints
//	for (int g=0; g<inst.numGen; g++) {
//		for (int t=0; t<nb_periods; t++) {
//			if (inst.must_run[g]) {
//				model.add( p_var[g][t] == capacity[g][t] - minimum_production_req[g] );
//			} else {
//				model.add( p_var[g][t] <= (capacity[g][t] - minimum_production_req[g]) * x[g][t] );
//			}
//		}
//	}
//
//	// production amounts
//	for (int g=0; g<inst.numGen; g++) {
//		for (int t=0; t<nb_periods; t++) {
//			model.add( p[g][t] == p_var[g][t] + x[g][t] * minimum_production_req[g] );
//		}
//	}
//
//	// ramp up constraints
//	for (int g=0; g<inst.numGen; g++) {
//		for (int t=1; t<nb_periods; t++) {
//			model.add( p_var[g][t] - p_var[g][t-1] <= inst.ramp_u_lim[g] * (double)period_len * x[g][t] - minimum_production_req[g] * s[g][t]);
//		}
//	}
//
//	// ramp down constraints
//	for (int g=0; g<inst.numGen; g++) {
//		for (int t=1; t<nb_periods; t++) {
//			model.add( p_var[g][t-1] - p_var[g][t] <= inst.ramp_d_lim[g] * (double)period_len * x[g][t-1] - minimum_production_req[g] * z[g][t] );
//		}
//	}
//
//	/** Objective Function **
//	IloExpr obj_func (env);
//
//	for (int g=0; g<inst.numGen; g++) {
//		for (int t=0; t<nb_periods; t++) {
//			obj_func += inst.start_cost[g] * s[g][t];								// start up cost
//			obj_func += minimum_production_req[g] * inst.var_cost[g] * x[g][t];		// cost of producing minimum production amount
//			obj_func += inst.no_load_cost[g] * (double)period_len / 60.0 * x[g][t];	// no-load cost
//			obj_func += inst.var_cost[g] * p_var[g][t];								// variable cost
//		}
//	}
//
//	// model-dependent variables, obj functions, and constraints
//	if (model_type == System) {
//		IloNumVarArray L (env, nb_periods, 0, IloInfinity, ILOFLOAT);	// load shedding
//
//		/** Constraints **
//		// aggregated demand constraints
//		for (int t=0; t<nb_periods; t++) {
//			IloExpr expr (env);
//			for (int g=0; g<inst.numGen; g++)
//				expr += p[g][t];
//			expr += L[t];
//			model.add( IloRange (env, aggregated_demand[t], expr) );
//		}
//
//		/** Objective Function **/
//		for (int t=0; t<nb_periods; t++) {
//			obj_func += inst.load_shedding_penalty * L[t];
//		}
//	}
//	else if (model_type == Transmission) {
//
//		/** Decision Variables **
//		IloArray< IloNumVarArray > L (env, inst.numBus);		// load shedding
//		IloArray< IloNumVarArray > T (env, inst.numBus);		// phase angles
//		for (int b=0; b<inst.numBus; b++) {
//			L[b] = IloNumVarArray(env, nb_periods, 0, IloInfinity, ILOFLOAT);
//			T[b] = IloNumVarArray(env, nb_periods, 0, 2*inst.pi, ILOFLOAT);
//
//			char buffer[30];
//
//			sprintf(buffer, "L_%d", b);
//			L[b].setNames(buffer);
//
//			sprintf(buffer, "T_%d", b);
//			T[b].setNames(buffer);
//		}
//
//		IloArray< IloNumVarArray > F (env, inst.numLine);		// flows
//		for (int l=0; l<inst.numLine; l++) {
//			F[l] = IloNumVarArray(env, nb_periods,
//					inst.min_arc_flow[l]*(double)period_len/60.0,
//					inst.max_arc_flow[l]*(double)period_len/60.0, ILOFLOAT);
//
//			char buffer[30];
//
//			sprintf(buffer, "F_%d", l);
//			F[l].setNames(buffer);
//		}
//
//		/** Constraints **/
//		// DC-approximation to AC flow
//		for (int l=0; l<inst.numLine; l++) {
//			for (int t=0; t<nb_periods; t++) {
//				int orig = inst.arcs[l].first;
//				int dest = inst.arcs[l].second;
//
//				model.add( F[l][t] == inst.susceptance[l] * (T[orig][t] - T[dest][t]) );
//			}
//		}
//
//		// Flow balance
//		for (int b=0; b<inst.numBus; b++) {
//			for (int t=0; t<nb_periods; t++) {
//
//				IloExpr expr (env);
//
//				// production
//				for (int g=0; g<inst.numGen; g++) {
//					if (inst.generator_loc[g] == b) expr += p[g][t];
//				}
//				// incoming arcs
//				for (int l=0; l<inst.numLine; l++) {
//					if (inst.arcs[l].second == b)   expr += F[l][t];
//				}
//				// outgoing arcs
//				for (int l=0; l<inst.numLine; l++) {
//					if (inst.arcs[l].first == b)    expr -= F[l][t];
//				}
//				// load shedding
//				expr += L[b][t];
//
//				// constraint
//				model.add( expr == demand[b][t] );
//
//				// clean expr from the memory
//				expr.end();
//			}
//		}
//
//		// remaining obj function components
//		for (int b=0; b<inst.numBus; b++) {
//			for (int t=0; t<nb_periods; t++) {
//				obj_func += inst.load_shedding_penalty * L[b][t];
//			}
//		}
//	}
//	else {
//		cout << "Invalid model number" << endl;
//	}
//
//	model.add( IloMinimize(env, obj_func) );
//
//	cplex.extract(model);
//	cplex.setParam(IloCplex::Threads, 1);
//	cplex.setParam(IloCplex::EpGap, 1e-2);
//}
//
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
