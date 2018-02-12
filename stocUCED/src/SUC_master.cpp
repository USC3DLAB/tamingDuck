//
//  SUC_master.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "SUC_master.hpp"

extern runType runParam;
extern ofstream optLog;

/****************************************************************************
 * Incumbent Callback
 * Used for retrieving expected initial-generation amounts for the optimal
 * master solution.
 ****************************************************************************/
IloCplex::Callback IncCallback(IloEnv env, SUCmaster &me) {
	return (IloCplex::Callback(new(env) SUCmaster::IncCallbackI(env, me)));
}

void SUCmaster::IncCallbackI::main() {
	
	// if the incumbent objective is improving,
	if ( getObjValue() < getIncumbentObjValue() ) {
		
		// set the generation amounts
		for (int g=0; g<me.numGen; g++) {
			me.setGenProd(g, 0, me.sub.expInitGen[g] );
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

void SUCmaster::LazySepCallbackI::main()
{
	// get the solution
	vector< vector<bool> > state;
	resize_matrix(state, me.numGen, me.numPeriods);

    IloNumArray x_vals (me.env);
	for (int g=0; g<me.numGen; g++) {
		getValues(x_vals, me.x[g]);
		for (int t=0; t<me.numPeriods; t++) {
			(x_vals[t] >= 0.5) ? state[g][t] = true : state[g][t] = false;
		}
	}
    x_vals.end();

	// set the master solution in the subproblems
	me.sub.setMasterSoln(state);
	
	// solve the subproblems
	bool isFeasible = me.sub.solve();
	
    // check if we need to add a cut (feasible and not within optimality tolerances, or infeasible)
    bool addCut = true;
    if (isFeasible) {
		double abs_diff = fabs(me.sub.getRecourseObjValue() - 1.0/((double)me.eta.getSize())*getValue(IloSum(me.eta)));
        if ( abs_diff/(fabs(me.sub.getRecourseObjValue())+1e-14) < me.cplex.getParam(IloCplex::EpGap)
            || abs_diff < me.cplex.getParam(IloCplex::EpAGap) ) {
            addCut = false;
        }
    }

	if (!addCut)	return;
	
	if (isFeasible) {
		/* optimality cut */
		
		if ( me.multicut ) {
			int optcut_cnt = 0;
			for (int s=0; s<me.eta.getSize(); s++) {
				
				double abs_diff = fabs(me.sub.objValues[s] - getValue(me.eta[s]));
				
				if ( abs_diff/(fabs(me.sub.getRecourseObjValue())+1e-14) < me.cplex.getParam(IloCplex::EpGap)
					|| abs_diff < me.cplex.getParam(IloCplex::EpAGap) ) {
					continue;
				}
				
				optcut_cnt++;
				
				// get the Benders's cut coefs
				SUCsubprob::BendersCutCoefs *cutCoefs = &(me.sub.multicutCoefs[s]);
				
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
				BendersCut = IloRange(me.env, cutCoefs->pi_b, me.eta[s]-pi_Tx);	// optimality cut
				
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
			
			optLog << "(optcut " << optcut_cnt << "/" << me.eta.getSize() << ")" << endl;
		}
		else
		{
			double sceProb = 1.0/(double)me.eta.getSize();
			
			IloExpr pi_Tx (me.env);
			double  pi_b = 0;
			
			for (int s=0; s<me.eta.getSize(); s++) {
				// get the Benders's cut coefs
				SUCsubprob::BendersCutCoefs *cutCoefs = &(me.sub.multicutCoefs[s]);
				
				// create the cut
				for (int g=0; g<me.numGen; g++) {
					for (int t=0; t<me.numPeriods; t++) {
						if ( fabs(cutCoefs->pi_T[g][t]) > 1e-10 ) {
							pi_Tx += cutCoefs->pi_T[g][t] * me.x[g][t];
						}
					}
				}
				pi_b += cutCoefs->pi_b;
			}
			
			pi_Tx *= sceProb;
			pi_b  *= sceProb;
			
			
			IloRange BendersCut;
			BendersCut = IloRange(me.env, pi_b, me.eta[0]-pi_Tx);	// optimality cut
			
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
	else {
		/* feasibility cut */
		optLog << "(feascut)" << endl;

		// get the Benders's cut coefs
		SUCsubprob::BendersCutCoefs *cutCoefs = &(me.sub.multicutCoefs[0]);
		
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
					expCapacity[g][t] += min( inst->observations["RT"].vals[rep][period][it->second], genPtr->maxCapacity );
				}
				else {
					it = inst->simulations.mapVarNamesToIndex.find(genPtr->name);
					for (int s=0; s<inst->simulations.vals.size(); s++) {
						expCapacity[g][t] += 1.0/(double)(inst->simulations.vals.size()) * min( inst->simulations.vals[s][period][it->second], genPtr->maxCapacity );
					}
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
			
			if ( t==0 && probType==ShortTerm ) {
				it = inst->observations["RT"].mapVarNamesToIndex.find( num2str(busPtr->regionId) );
				busLoad[b][t] = inst->observations["RT"].vals[rep][period][it->second] * busPtr->loadPercentage;
			}
			else {
				it = inst->observations["DA"].mapVarNamesToIndex.find( num2str(busPtr->regionId) );
				busLoad[b][t] = inst->observations["DA"].vals[rep][period][it->second] * busPtr->loadPercentage;
			}
		}
	}
	
	/* Spinning Reserve */
	for (int b=0; b<numBus; b++) {
		for (int t=0; t<numPeriods; t++) {
			busLoad[b][t] *= (1+spinReservePerc);
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
	multicut		= true;
	
	/* Prepare Model-Dependent Input Data */
	preprocessing();

	/* Master Formulation */
	// create the variables
	//IloNumVar L (env);
	
	eta = IloNumVarArray (env, multicut ? runParam.numLSScen : 1, 0, IloInfinity, ILOFLOAT);	// 2nd-stage approximation
	
	IloArray<IloNumVarArray> p (env, numGen);	// production amounts
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
    
    /* demand-based valid inequality
    // the capacities of operational generators must exceed the system demand at any point
    for (int t=0; t<numPeriods; t++) {
        IloExpr expr (env);
        for (int g=0; g<numGen; g++) {
            expr += expCapacity[g][t] * x[g][t];
        }
		expr += L;
        model.add( expr >= sysLoad[t] );
        expr.end();
    }
	*/
	
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
		
		int t=0;
		if ( beginMin != 0 ) {
			IloConstraint c ( p[g][t] - getEDGenProd(g, t-1) <= genPtr->rampUpLim*periodLength );
			sprintf(buffer, "RU_%d_%d", g, t);
			c.setName(buffer);
			model.add(c);
		}
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
		
		int shutDownPeriod = -9999;
		if ( probType == ShortTerm && genPtr->isDAUCGen && beginMin != 0 ) {
			// check if this generator was producing at t-1
			bool wasProducing = (getEDGenProd(g, -1) > 1e-8);
			
			if (wasProducing) {
				// check if this generator is shutting down at some point in ST-UC planning horizon
				bool willShutDown = false;
				int t;
				for (t=0; t<numPeriods; t++) {
					if (!getGenState(g, t)) {
						willShutDown = true;
						break;
					}
				}
				
				if (willShutDown) {
					// check if it can ramp down
					if ( genPtr->rampDownLim*periodLength*(double)(t+1) < getEDGenProd(g, -1) ) {
						shutDownPeriod = t;
					}
				}
			}
		}
		
		if (shutDownPeriod >= 0) {
			optLog << "Warning: Generator " << g << " (" << genPtr->name << ") cannot ramp down to 0 in the ST-UC problem" << endl;
		}
		
		int t=0;
		if ( beginMin != 0 ) {
			double rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			IloConstraint c ( getEDGenProd(g, t-1) - p[g][t] <= rampDownRate );
			sprintf(buffer, "RD_%d_%d", g, t); c.setName(buffer); model.add(c);
		}
		for (t=1; t<numPeriods; t++) {
			double rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			IloConstraint c ( p_var[g][t-1] - p_var[g][t] <= rampDownRate * x[g][t-1] - minGenerationReq[g] * z[g][t] );
			sprintf(buffer, "RD_%d_%d", g, t); c.setName(buffer); model.add(c);
		}
	}

	
	// create the objective function
	IloExpr obj (env);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			obj += genPtr->startupCost * s[g][t];					// start up cost
			obj += genPtr->noLoadCost*periodLength/60.0 * x[g][t];	// no-load cost
		}
	}
	for (int s=0; s<eta.getSize(); s++) {
		obj += 1.0/(double)eta.getSize() * eta[s];
	}
	//obj += loadShedPenaltyCoef * L;
	
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
		IloArray< IloNumVarArray > L (env, numBus);	// load shedding
		IloArray< IloNumVarArray > O (env, numBus);	// over generation
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
	
    // set the objective function
	model.add( IloMinimize(env, obj) );
	
    // formulate the subproblem
	sub.formulate(inst, probType, modelType, beginMin, rep);
	
	// formulate the warm-up problem
	warmUpProb.formulate(inst, probType, modelType, beginMin, rep);

    /*************************************************************************
     * DISABLED: Didn't improve the progress of the algorithm. Donno why..
     *
     * // compute a lower bound on eta
     * double eta_lb = sub.computeLowerBound();
     * eta.setLB( eta_lb );
     ************************************************************************/
     
    // prepare the solver
	cplex.extract(model);

    /*************************************************************************
     * DISABLED: Didn't improve the progress much.. reduced it in some cases.
     *
     * // assign priorities to state variables
     * for (int g=0; g<inst.G; g++) {
	 *    for (int t=0; t<inst.T; t++) {
     *         cplex.setPriority(x[g][t], 1);
     *     }
     * }
     ************************************************************************/
    
	cplex.use(LazySepCallback(env, *this));
	if (probType==DayAhead && beginMin == 0) {
		cplex.use(IncCallback(env, *this));
	}
	cplex.setOut(optLog);
	cplex.setWarning(optLog);
	cplex.setParam(IloCplex::Threads, 1);
//	cplex.setParam(IloCplex::MIPEmphasis, IloCplex::MIPEmphasisFeasibility);
//    cplex.setParam(IloCplex::FPHeur, 2);
//    cplex.setParam(IloCplex::HeurFreq, 10);
//    cplex.setParam(IloCplex::LBHeur, 1);    // ??
    cplex.setParam(IloCplex::EpGap, 1e-2);
//    cplex.setParam(IloCplex::Reduce, 0);    // due to callbacks
//    cplex.setParam(IloCplex::PreLinear, 0); // due to callbacks
	cplex.setParam(IloCplex::TiLim, (probType == DayAhead)*7200 + (probType == ShortTerm)*600);
//	cplex.setParam(IloCplex::TiLim, (probType == DayAhead)*150 + (probType == ShortTerm)*60);
}

bool SUCmaster::solve () {
	try{
		bool status;
		
		/* warm up *
		optLog << "Executing warm up MIP..." << endl;
		status = warmUpProb.solve();
		if (status) {
			setWarmUp();
			optLog << "Warm up model provided " << warmUpProb.cplex.getSolnPoolNsolns() << " solutions (Best Obj= " << warmUpProb.getObjValue() << ")." << endl;
		} else {
			optLog << "Warm up has failed." << endl;
		}
		/*        */
		
		// Benders' decomposition
		optLog << "Executing Benders' decomposition..." << endl;
		status = cplex.solve();
		if (status) {
			optLog << "Optimization is completed with status " << cplex.getCplexStatus() << endl;
			optLog << "Obj = \t" << cplex.getObjValue() << endl;
			optLog << "LB = \t" << cplex.getBestObjValue() << endl;
		} else {
			cplex.exportModel("infeasible.lp");
			optLog << "Benders' decomposition has failed." << endl;
		}
		
		// Record the optimal solution
		if (status) {
			for (int g=0; g<numGen; g++) {
				Generator *genPtr = &(inst->powSys->generators[g]);
				
				for (int t=0; t<numPeriods; t++) {
					if ( (probType == DayAhead && genPtr->isDAUCGen) || (probType == ShortTerm && !genPtr->isDAUCGen) ) {
						setGenState(g,t, cplex.getValue(x[g][t]));
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
		cplex.addMIPStart(vars, vals);
		
		vars.end();
		vals.end();
	}
}

/*
void master::warmStart(Solution &soln)
{
    this->soln = soln;
 
    IloNumVarArray vars (env);
    IloNumArray vals (env);
 
    for (int g=0; g<inst->G; g++) {
        for (int t=0; t<inst->T; t++) {
            vars.add( s[g][t] );
            vals.add( soln.s[g][t] );
            
            vars.add( x[g][t] );
            vals.add( soln.x[g][t] );

            vars.add( z[g][t] );
            vals.add( soln.z[g][t] );
        }
    }
    
    cplex.addMIPStart(vars, vals);
    
    vars.end();
    vals.end();
}

void master::warmStart(vector<Solution> &solnPool)
{
    this->soln = solnPool[0];
    
    IloNumVarArray vars (env);
    IloNumArray vals (env);
    
    for (int sol=0; sol<solnPool.size(); sol++)
    {
        for (int g=0; g<inst->G; g++) {
            for (int t=0; t<inst->T; t++) {
                vars.add( s[g][t] );
                vals.add( solnPool[sol].s[g][t] );
                
                vars.add( x[g][t] );
                vals.add( solnPool[sol].x[g][t] );
                
                vars.add( z[g][t] );
                vals.add( solnPool[sol].z[g][t] );
            }
        }
        cplex.addMIPStart(vars, vals);
    }
    
    vars.end();
    vals.end();
}

*/

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
	else if (reqSolnComp < inst->solution.x[genId].size()) {	// return the corresponding solution
		return round(inst->solution.x[genId][reqSolnComp]);
	}
	else {														// asking what's beyond the planning horizon, we return the last solution
		return round(inst->solution.x[genId][ inst->solution.x[genId].size()-1 ]);
	}
}

/****************************************************************************
 * setGenState
 * - Fills the (genId, correspondingComponent) of the Solution.x object.
 ****************************************************************************/
void SUCmaster::setGenState(int genId, int period, double value) {
	// which Solution component is being set?
	int solnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// set the solution
	if (solnComp >= 0 && solnComp < inst->solution.x[genId].size()) {
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
	else if (reqSolnComp < inst->solution.x[genId].size()) {	// return the corresponding solution
		return inst->solution.g_ED[genId][reqSolnComp];
	}
	else {														// asking what's beyond the planning horizon, we return the last solution
		cout << "Error: Production levels beyond the planning horizon are not available" << endl;
		exit(1);
	}
}

/****************************************************************************
 * setGenProd
 * - Fills the (genId, correspondingComponent) of the Solution.g object.
 ****************************************************************************/
void SUCmaster::setGenProd(int genId, int period, double value) {
	// which Solution component is being set?
	int solnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// set the solution
	if (solnComp >= 0 && solnComp < inst->solution.g_UC[genId].size()) {
		for (int t=solnComp; t<solnComp+numBaseTimePerPeriod; t++) {
			inst->solution.g_UC[genId][t] = value;
		}
	}
	else {
		// Setting generator production at time that is beyond the planning horizon
	}
}
