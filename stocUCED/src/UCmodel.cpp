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
 * - Initializes basic parameters
 * - Fills the following containers, according to model units & assumptions:
 *	 minGenerationReq
 *   minUpTimePeriods
 *   minDownTimePeriods
 *   capacity
 *   busLoads
 *   aggLoads
 * - If it is a DA-UC problem, sets the capacity and minGenerationReq to 0
 * for generators which must be scheduled at ST-UC.
 *
 * - Assumption 1: Min generation requirement must be less than ramping
 * rates, otherwise, generators cannot switch on/off.
 * - Assumption 2: Min up/down times must at least be 1 period.
 ****************************************************************************/
void UCmodel::preprocessing ()
{
	/* basic parameters */
	numGen = inst->powSys->numGen;
	numLine = inst->powSys->numLine;
	numBus = inst->powSys->numBus;
	numBatteries = inst->powSys->numBatteries;
	
	numPeriods	 = (probType == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
	periodLength = (probType == DayAhead) ? runParam.DA_resolution : runParam.ST_resolution;
	
	numBaseTimePerPeriod = (int)round(periodLength) / runParam.ED_resolution;
	
	
	/* initialize containers */
	minGenerationReq.resize(numGen);					// minimum production requirements
	minUpTimePeriods.resize(numGen);					// minimum uptime in periods
	minDownTimePeriods.resize(numGen);					// minimum downtime in periods
	resize_matrix(capacity, numGen, numPeriods);		// generator capacities
	sysLoad.resize(numPeriods);							// aggregated system load
	resize_matrix(busLoad, numBus, numPeriods);			// individual bus loads
	
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
	
	/* Mean Generator Capacity */
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);

		auto it = inst->actuals.mapVarNamesToIndex.find(genPtr->name);
		if ( it != inst->actuals.mapVarNamesToIndex.end() ) {
			/* supply info found within the time series */
			
			int period;
			for (int t=0; t<numPeriods; t++) {
				period = (beginMin/periodLength)+(t*numBaseTimePerPeriod);

				if (t == 0 && probType == ShortTerm) {
					// use actual values for the 1st-period of the ST-UC problem
					it = inst->actuals.mapVarNamesToIndex.find(genPtr->name);
					capacity[g][t] = min(inst->actuals.vals[rep][period][it->second], genPtr->maxCapacity);
				}
				else {
					if (runParam.updateForecasts && probType != DayAhead) {
						// use updated forecasts (only) for the ST-UC problem
						it = inst->meanForecast["RT"].mapVarNamesToIndex.find(genPtr->name);
						capacity[g][t] = min(inst->meanForecast["RT"].vals[rep][period][it->second], genPtr->maxCapacity);
					}
					else {
						// use DA forecasts if no updates are made, or for the DA-UC problem
						it = inst->meanForecast["DA"].mapVarNamesToIndex.find(genPtr->name);
						capacity[g][t] = min(inst->meanForecast["DA"].vals[rep][period][it->second], genPtr->maxCapacity);
					}
				}
			}
		}
		else {
			/* deterministic supply */
			fill(capacity[g].begin(), capacity[g].end(), genPtr->maxCapacity);
		}
	}

	/* Misc: Generators which are scheduled in ST-UC, must have 0 min-generation-requirement & capacity in DA-UC problem */
	if (probType == DayAhead) {
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst->powSys->generators[g]);
			
			if (!genPtr->isDAUCGen) {
				minGenerationReq[g] = 0.0;
				fill(capacity[g].begin(), capacity[g].end(), 0.0);
			}
		}
	}
	
	/* Load */
	for (int b=0; b<numBus; b++) {
		Bus *busPtr = &(inst->powSys->buses[b]);
		
		auto it = inst->actuals.mapVarNamesToIndex.find( num2str(busPtr->regionId) );
		
		int period;
		for (int t=0; t<numPeriods; t++) {
			period = (beginMin/periodLength)+(t*numBaseTimePerPeriod);
			
			if ( probType == DayAhead ) {
				// use DA forecast of the load in the DA-UC problem
				it = inst->meanForecast["DA"].mapVarNamesToIndex.find( num2str(busPtr->regionId) );
				busLoad[b][t] = inst->meanForecast["DA"].vals[rep][period][it->second] * busPtr->loadPercentage;
			}
			else {
				// use actual values of the load in the ST-UC problems (assuming the forecasts are accurate, which is quite true)
				it = inst->actuals.mapVarNamesToIndex.find( num2str(busPtr->regionId) );
				busLoad[b][t] = inst->actuals.vals[rep][period][it->second] * busPtr->loadPercentage;
			}
		}
	}

	/* Spinning Reserve */
	for (int b=0; b<numBus; b++) {
		for (int t=0; t<numPeriods; t++) {
			busLoad[b][t] *= (1+runParam.resPerc_UC);
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

void UCmodel::initializeVariables () {
	// Generator dvars
	s	  = IloArray<IloNumVarArray> (env, numGen);	// start up vars
	x	  = IloArray<IloNumVarArray> (env, numGen);	// state vars
	z	  = IloArray<IloNumVarArray> (env, numGen);	// shut down vars
	p	  = IloArray<IloNumVarArray> (env, numGen);	// production amounts
	p_var = IloArray<IloNumVarArray> (env, numGen);	// variable-production amounts	
	// Battery dvars
	v_pos = IloArray<IloNumVarArray> (env, numBatteries);	// charge amounts
	v_neg = IloArray<IloNumVarArray> (env, numBatteries);	// discharge amounts
	I     = IloArray<IloNumVarArray> (env, numBatteries);	// battery state

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
	
	for (int bt=0; bt<numBatteries; bt++) {
		v_pos[bt] = IloNumVarArray(env, numPeriods, 0, inst->powSys->batteries[bt].maxCapacity, ILOFLOAT);
		v_neg[bt] = IloNumVarArray(env, numPeriods, 0, inst->powSys->batteries[bt].maxCapacity, ILOFLOAT);
		I[bt] = IloNumVarArray(env, numPeriods, 0, inst->powSys->batteries[bt].maxCapacity, ILOFLOAT);
		
		sprintf(buffer, "vp_%d", bt);
		v_pos[bt].setNames(buffer);
		sprintf(buffer, "vn_%d", bt);
		v_neg[bt].setNames(buffer);
		sprintf(buffer, "I_%d", bt);
		I[bt].setNames(buffer);
		
		model.add(v_pos[bt]);
		model.add(v_neg[bt]);
		model.add(I[bt]);
	}
}

/****************************************************************************
 * formulate
 * - Performs preprocessing on the instance, to compute necessary model
 * parameters.
 * - Formulates the mathematical program.
 ****************************************************************************/
void UCmodel::formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep) {
	
	/* Initializations */
	this->inst		= &inst;
	this->beginMin	= beginMin;
	this->probType	= probType;
	this->modelType = modelType;
	this->rep		= rep;

	/* Prepare Model-Dependent Input Data */
	preprocessing();
	
	/* Prepare the Mathematical Model */
	initializeVariables();
	
	/**** Constraints (Traditional Formulation) ****/
	
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
	
	// minimum uptime constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		if ( (probType == DayAhead && !genPtr->isDAUCGen) || (probType == ShortTerm && genPtr->isDAUCGen) ) {
			continue;	// skip scheduling constraints for not-to-be-scheduled generators
		}

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
	
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			if (genPtr->isMustUse) {
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
		
		// ignore ramping constraints for generators that are not meant to be scheduled in DA-UC problem
		if (probType == DayAhead && !genPtr->isDAUCGen)	continue;
		
		/* t = 0 */
		int t = 0;
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
		
		// ignore ramping constraints for generators that are not meant to be scheduled in DA-UC problem
		if (probType == DayAhead && !genPtr->isDAUCGen)	continue;

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
	
	// battery state
	for (int bt=0; bt<numBatteries; bt++) {
		Battery* bat_ptr = &inst.powSys->batteries[bt];
		
		double dissipationCoef = pow(bat_ptr->dissipationCoef, periodLength/60.0);
		double chargingLossCoef = pow(bat_ptr->chargingLossCoef, periodLength/60.0);
		double dischargingLossCoef = pow(bat_ptr->dischargingLossCoef, periodLength/60.0);

		int t=0;
		IloConstraint c ( I[bt][t] == getBatteryState(bat_ptr->id, -1) * dissipationCoef + v_pos[bt][t] * chargingLossCoef - v_neg[bt][t] * dischargingLossCoef);
		sprintf(buffer, "Bt_%d_%d", bt, t); c.setName(buffer); model.add(c);

		for (t=1; t<numPeriods; t++) {
			IloConstraint c ( I[bt][t] == I[bt][t-1] * dissipationCoef + v_pos[bt][t] * chargingLossCoef - v_neg[bt][t] * dischargingLossCoef );
			sprintf(buffer, "Bt_%d_%d", bt, t); c.setName(buffer); model.add(c);
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

	
	/* ST-UC rampability constraints *
	if (probType == ShortTerm) {
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst.powSys->generators[g]);
			
			if (genPtr->type != Generator::SOLAR && genPtr->type != Generator::WIND && genPtr->isDAUCGen) {
				int t = numPeriods-1;
				
				int target_idx = 0;
				if (beginMin / runParam.ED_resolution + numPeriods * numBaseTimePerPeriod >= runParam.numPeriods ) {
					target_idx = numPeriods-1;
				} else {
					target_idx = numPeriods;
				}
				
				double target = getUCGenProd(g, target_idx);
				
				/* ramp-up */
				//				model.add( target - p[g][t] - delta_pos[g][t] <= genPtr->rampUpLim * runParam.DA_resolution );
				
				/* ramp-down */
				//				model.add( p[g][t] - target - delta_neg[g][t] <= genPtr->rampDownLim * runParam.DA_resolution );
//			}
//		}
//	}
	/* */

	
	/** Objective Function **/
	IloExpr obj (env);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		double varCostPerPeriod    = genPtr->variableCost * periodLength/60.0;
		double noloadCostPerPeriod = genPtr->noLoadCost * periodLength/60.0;
		
		for (int t=0; t<numPeriods; t++) {
			obj += genPtr->startupCost * s[g][t];						// start up cost
			obj += minGenerationReq[g] * varCostPerPeriod * x[g][t];	// cost of producing minimum generation requirements
			obj += varCostPerPeriod * p_var[g][t];						// variable cost
			obj += noloadCostPerPeriod * x[g][t];						// no-load cost
		}
	}	
	
	/** Model-dependent obj function components, constraints **/
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
			for (int bt=0; bt<numBatteries; bt++) {
				expr -= v_pos[bt][t] + v_neg[bt][t];
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
				for (int g=0; g < (int) busPtr->connectedGenerators.size(); g++) {
					expr += p[ busPtr->connectedGenerators[g]->id ][t];
				}
				
				// storage (iterate over connected batteries)
				for (int bt=0; bt < (int) busPtr->connectedBatteries.size(); bt++) {
					expr -= v_pos[ busPtr->connectedBatteries[bt]->id ][t] + v_neg[ busPtr->connectedBatteries[bt]->id ][t];
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
		
		// load-shedding penalties
		for (int b=0; b<numBus; b++) {
			for (int t=0; t<numPeriods; t++) {
				obj += loadShedPenaltyCoef * L[b][t];
				obj += overGenPenaltyCoef * O[b][t];
			}
		}
	}
	
	/** Finalize **/
	model.add( IloMinimize(env, obj) );
	obj.end();
	
	cplex.extract(model);
	cplex.setParam(IloCplex::EpGap, 1e-2);
	cplex.setOut(inst.out());
	cplex.setWarning(inst.out());
}

void UCmodel::saveSolution() {
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			if ( (probType == DayAhead && genPtr->isDAUCGen) || (probType == ShortTerm && !genPtr->isDAUCGen) ) {
				setGenState(g,t, cplex.getValue(x[g][t]));
			}
			// Note: Production = 0 when generator is not operational. Below line prevents numerical errors
			// where there is >0 production but =0 commitment.
			setUCGenProd(g, t, cplex.getValue(p[g][t]) * getGenState(g, t));
		}
	}
	
	double totLoadShed = 0;
	for (int b=0; b<numBus; b++) {
		for (int t=0; t<numPeriods; t++) {
			if ( cplex.getValue(L[b][t]) > 1e-6 ) {
				totLoadShed += cplex.getValue(L[b][t]);
			}
		}
	}
	if (totLoadShed > 0) {
		cout << "[LS! " << setprecision(1) << totLoadShed << " MWs] ";
	}
}

/****************************************************************************
 * solve
 * - Solves the model and saves the solution into inst->solution depending
 * on the flag.
 ****************************************************************************/
bool UCmodel::solve(bool saveSol) {
	bool status;
	try {
		status = cplex.solve();
		
		// record the solution
		if (status && saveSol) {
			saveSolution();
		}
		
		if (!status) {
			cplex.exportModel("infeasible_UC.lp");
			exit(5);
		}
	}
	catch (IloException &e) {
		cout << e << endl;
	}
	
	return status;
}

/****************************************************************************
 * solve
 * - Solves the model and saves the solution into inst->solution.
 ****************************************************************************/
bool UCmodel::solve() {
	return solve(true);
}

/****************************************************************************
 * getGenState
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded status of the generator.
 * - If available, the function returns info prior to planning horizon. If
 * info is not available, the function quits.
 ****************************************************************************/
bool UCmodel::getGenState(int genId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// return the requested generator state
	if (reqSolnComp < 0) {
		if (runParam.useGenHistory) {
			// find the requested day's solution
			list<Solution>::reverse_iterator rit = inst->solList.rbegin();
			while (rit != inst->solList.rend()) {
				reqSolnComp += runParam.numPeriods;
				if (reqSolnComp >= 0) break;
				
				++rit;
			}
			
			// if found, return, otherwise return true.
			if (rit != inst->solList.rend()) {
				return rit->x[genId][reqSolnComp];
			}
			else {
				return true;	// if not found, assume that the generator has been operational for a long time
			}
		}
		else {
			return true;	// assume the generator has been operational for a long time
		}
	}
	else if (reqSolnComp < (int) inst->solution.x[genId].size()) {	// return the corresponding solution
		return ( inst->solution.x[genId][reqSolnComp] > EPSzero );
	}
	else {														// asking what's beyond the planning horizon, we return the last solution
		return ( inst->solution.x[genId][ inst->solution.x[genId].size()-1 ] > EPSzero );
	}
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
int UCmodel::checkShutDownRampDownInconsistency (int g) {
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
 * getGenProd
 * - If past info is requested, we will return past ED generation info if it
 * is available AND we are allowed to use it (due to configurations).
 * - If past info is requested, but either we don't have it or not allowed
 * to use it, we use the 1st-period DA-UC generations of DA generators as
 * history in the ST-UC problem.
 * - In all other times, we return ED production levels.
 ****************************************************************************/
double UCmodel::getGenProd(int g, int t) {
	if (t < 0 && beginMin == 0) {									// info requested is prior to the planning horizon
		if ( runParam.useGenHistory && inst->solList.size() > 0 ) {	// use history
			return getEDGenProd(g, -1);	// this will return the final ED gen levels from the prev day sol
		}
		else if (probType == ShortTerm && inst->powSys->generators[g].isDAUCGen){
			return getUCGenProd(g, 0);
		}
	} else {
		return getEDGenProd(g, t);
	}
	return -INFINITY;
}

double UCmodel::getBatteryState(int batteryId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// return the requested generator state
	if (reqSolnComp > 0) {
		return inst->solution.btState_ED[batteryId][reqSolnComp];
	} else {
		if (runParam.useGenHistory && inst->solList.size() > 0) {
			return inst->solList.back().btState_ED[batteryId][runParam.numPeriods-1];	// the latest battery state
		} else {
			return 0.0;
		}
	}
}

/****************************************************************************
 * getUCGenProd
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded generation of the generator by the DA model.
 ****************************************************************************/
double UCmodel::getUCGenProd(int genId, int period) {
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

/****************************************************************************
 * getEDGenProd
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded generation of the generator by the ED model.
 * - If past info is requested, it will return the final ED generation of the
 * previous day.
 ****************************************************************************/
double UCmodel::getEDGenProd(int genId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// return the requested generator state
	if (reqSolnComp < 0) {										// all generators are assumed to be online, for a long time, at t=0.
		if (inst->solList.size() > 0) {
			return inst->solList.back().g_ED[genId][runParam.numPeriods-1];
		}
		else {
			cout << "Error: Initial production levels are not available" << endl;
			exit(1);
		}
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
 * setGenState
 * - Fills the (genId, correspondingComponent) of the Solution.x object.
 ****************************************************************************/
void UCmodel::setGenState(int genId, int period, double value) {
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
 * setUCGenProd
 * - Fills the (genId, correspondingComponent) of the Solution.g object.
 ****************************************************************************/
void UCmodel::setUCGenProd(int genId, int period, double value) {
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
 * getObjValue
 * - Returns the objective value of the last optimization.
 ****************************************************************************/
double UCmodel::getObjValue() {
	return cplex.getObjValue();
}
