//
//  subprob.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "SUC_subprob.hpp"

extern runType runParam;

SUCsubprob::SUCsubprob () {
	model = IloModel (env);
	cplex = IloCplex (env);

	cons  = IloRangeArray (env);
	duals = IloNumArray (env);

	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::PreInd, 0);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);	// Important
}

SUCsubprob::~SUCsubprob () {
	env.end();
}

/****************************************************************************
 * preprocessing
 * - Initializes basic parameters
 * - Fills the following containers, according to model units & assumptions:
 *	 minGenerationReq
 *   minDownTimePeriods
 *   expCapacity
 *   busLoads
 *   aggLoads
 * - If it is a DA-UC problem, sets the capacity and minGenerationReq to 0
 * for generators which must be scheduled at ST-UC.
 *
 * - Assumption 1: Min generation requirement must be less than ramping
 * rates, otherwise, generators cannot switch on/off.
 * - Assumption 2: Min up/down times must at least be 1 period.
 ****************************************************************************/
void SUCsubprob::preprocessing ()
{
	/* basic parameters */
	numGen	   = inst->powSys->numGen;
	numLine    = inst->powSys->numLine;
	numBus     = inst->powSys->numBus;
	
	numPeriods	 = (probType == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
	periodLength = (probType == DayAhead) ? runParam.DA_resolution : runParam.ST_resolution;
	
	numBaseTimePerPeriod = (int)round(periodLength) / runParam.ED_resolution;
	
	/* initialize containers */
	cutCoefs.initialize(numGen, numPeriods);
	minGenerationReq.resize(numGen);					// minimum production requirements
	maxCapacity.resize(numGen);							// maximum generator capacities
	if (modelType == System) { sysLoad.resize(numPeriods); }					// aggregated system load
	else					 { resize_matrix(busLoad, numBus, numPeriods); }	// individual bus loads
	
	/* Min Generation Amounts */
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		minGenerationReq[g] = genPtr->minGenerationReq;
		if (minGenerationReq[g] > min(genPtr->rampUpLim * periodLength, genPtr->rampDownLim * periodLength)) {
			minGenerationReq[g] = min(genPtr->rampUpLim * periodLength, genPtr->rampDownLim * periodLength);
		}
	}

	/* Max Generator Capacity */
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		maxCapacity[g] = genPtr->maxCapacity;
	}
	
	/* Misc: Generators which are scheduled in ST-UC, must have 0 min-generation-requirement & capacity in DA-UC problem */
	if (probType == DayAhead) {
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst->powSys->generators[g]);
			
			if (!genPtr->isDAUCGen) {
				minGenerationReq[g] = 0.0;
				maxCapacity[g] = 0.0;
			}
		}
	}
	
	/* Load */
	auto dataPtr = (probType == DayAhead) ? &(inst->detObserv[0]) : &(inst->detObserv[1]);
	if (modelType == System) {
		fill( sysLoad.begin(), sysLoad.end(), 0.0 );		// initialize to 0
		for (int t=0; t<numPeriods; t++) {
			for (int r=0; r<dataPtr->numVars; r++) {
				sysLoad[t] += dataPtr->vals[0][t*numBaseTimePerPeriod][r];
			}
		}
	}
	else {
		for (int b=0; b<numBus; b++) {
			Bus *busPtr = &(inst->powSys->buses[b]);
			int r = dataPtr->mapVarNamesToIndex[ num2str(busPtr->regionId) ];
			
			for (int t=0; t<numPeriods; t++) {
				busLoad[b][t] = dataPtr->vals[0][t*numBaseTimePerPeriod][r] * busPtr->loadPercentage;
			}
		}
	}
	
	/* Uncertainty related */
	scenSet = (probType == DayAhead) ? &(inst->stocObserv[0]) : &(inst->stocObserv[1]);
	numScen = (int)scenSet->vals.size();
	
	sceProb.resize(numScen);
	fill(sceProb.begin(), sceProb.end(), 1.0/(double)numScen);		// equal probability scenarios
	
	// TODO: remove this
	cout << "Solving for " << numScen << " scenarios " << endl;
}

void SUCsubprob::formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin)
{
	// get a handle on the instance and model type
	this->inst		= &inst;
	this->beginMin	= beginMin;
	this->probType	= probType;
	this->modelType = modelType;
	
	preprocessing();
		
	// formulate the desired problem
	if (modelType == System) {
		formulate_aggregate_system();
	} else {
		formulate_nodebased_system();
	}
}

void SUCsubprob::formulate_production()
{
	/** Decision Variables **/
	p = IloArray< IloNumVarArray > (env, numGen);	// production amounts
	for (int g=0; g<numGen; g++) {
		p[g]	 = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
		
		sprintf(buffer, "p_%d", g);
		p[g].setNames(buffer);
	}

	/** Constraints **/
 
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			IloRange con;
			if (genPtr->isMustUse) {
                con = IloRange (env, genPtr->maxCapacity, p[g][t], genPtr->maxCapacity);
			} else {
                con = IloRange (env, -IloInfinity, p[g][t], genPtr->maxCapacity);
			}
			cons.add( con );
		}
	}
	
	// minimum production constraints
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++) {
			IloRange con (env, minGenerationReq[g], p[g][t]);
			cons.add( con );
		}
	}
	
	// ramp up constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++) {
			IloRange con (env, -IloInfinity, p[g][t] - p[g][t-1], genPtr->rampUpLim * periodLength);
			cons.add( con );
		}
	}
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++) {
			IloRange con (env, -IloInfinity, p[g][t-1] - p[g][t], genPtr->rampDownLim * periodLength);
			cons.add( con );
		}
	}
}

void SUCsubprob::formulate_nodebased_system()
{
	formulate_production();
	
	/** Decision Variables **/
    L = IloArray< IloNumVarArray > (env, numBus);	// load shedding
	IloArray< IloNumVarArray > O (env, numBus);		// over generation
	IloArray< IloNumVarArray > T (env, numBus);		// phase angles
	for (int b=0; b<numBus; b++) {
		L[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
		O[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
		T[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
		
		sprintf(buffer, "L_%d", b);
		L[b].setNames(buffer);
		
		sprintf(buffer, "O_%d", b);
		O[b].setNames(buffer);
		
		sprintf(buffer, "T_%d", b);
		T[b].setNames(buffer);
	}
	
	IloArray< IloNumVarArray > F (env, numLine);		// flows
	for (int l=0; l<numLine; l++) {
		F[l] = IloNumVarArray(env, numPeriods, -IloInfinity, IloInfinity, ILOFLOAT);
		
		sprintf(buffer, "F_%d", l);
		F[l].setNames(buffer);
	}
	
	/** Constraints **/
	// Flow upper and lower bounds
	for (int l=0; l<numLine; l++) {			// upper bound
		for (int t=0; t<numPeriods; t++) {
			IloRange con(env, -IloInfinity, F[l][t], inst->powSys->lines[l].maxFlowLim);
			cons.add(con);
		}
	}
	for (int l=0; l<numLine; l++) {			// lower bound
		for (int t=0; t<numPeriods; t++) {
			IloRange con(env, inst->powSys->lines[l].minFlowLim, F[l][t]);
			cons.add(con);
		}
	}
	
	// Phase-angle upper bounds
	for (int b=0; b<numBus; b++) {
		Bus *busPtr = &(inst->powSys->buses[b]);
		
		for (int t=0; t<numPeriods; t++) {
			IloRange con(env, -IloInfinity, T[b][t], busPtr->maxPhaseAngle - busPtr->minPhaseAngle);
			cons.add(con);
		}
	}
	
	// DC-approximation to AC flow
	for (int l=0; l<numLine; l++) {
		Line *linePtr = &(inst->powSys->lines[l]);
		
		int orig = linePtr->orig->id;
		int dest = linePtr->dest->id;

		for (int t=0; t<numPeriods; t++) {
			model.add( F[l][t] == linePtr->susceptance * (T[orig][t] - T[dest][t]) );
		}
	}

	// Flow balance
	for (int b=0; b<numBus; b++) {
		Bus *busPtr = &(inst->powSys->buses[b]);
		
		for (int t=0; t<numPeriods; t++) {
			IloExpr expr (env);
			
			// production (iterate over connected generators)
			for (int g=0; g<busPtr->connectedGenerators.size(); g++) {
				expr += p[ busPtr->connectedGenerators[g]->id ][t];
			}


			// in/out flows (iterate over all arcs)
			for (auto linePtr = inst->powSys->lines.begin(); linePtr != inst->powSys->lines.end(); ++linePtr) {
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
			IloRange con(env, busLoad[b][t], expr, busLoad[b][t]);
			cons.add(con);
			
			expr.end();
		}
	}
	
	/** Objective Function **/
	IloExpr obj (env);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			obj += genPtr->variableCost * periodLength/60.0 * p[g][t];
		}
	}
	
	for (int b=0; b<numBus; b++) {
		for (int t=0; t<numPeriods; t++) {
			obj += loadShedPenaltyCoef * L[b][t];
			obj -= overGenPenaltyCoef * O[b][t];
		}
	}
	
	/** Finalize **/
	model.add( IloMinimize(env, obj) );
	model.add( cons );
	cplex.extract(model);
}

void SUCsubprob::formulate_aggregate_system()
{
	formulate_production();

	/** Decision Variables **/
	IloNumVarArray L (env, numPeriods, 0, IloInfinity, ILOFLOAT);	// load shedding
	IloNumVarArray O (env, numPeriods, 0, IloInfinity, ILOFLOAT);	// over generation

	/** Constraints **/
	// aggregated demand constraints
	for (int t=0; t<numPeriods; t++) {
		IloExpr expr (env);
		for (int g=0; g<numGen; g++)	expr += p[g][t];
		expr += L[t];
		expr -= O[t];
		IloRange con (env, sysLoad[t], expr, sysLoad[t]);
		cons.add( con );
	}
	
	/** Objective Function **/
	IloExpr obj (env);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);

		for (int t=0; t<numPeriods; t++) {
			obj += genPtr->variableCost * periodLength/60.0 * p[g][t];
		}
	}
	
	for (int t=0; t<numPeriods; t++) {
		obj += loadShedPenaltyCoef * L[t];
		obj -= overGenPenaltyCoef * O[t];
	}
	
	// prepare the model and the solver
	model.add( IloMinimize(env, obj) );
	model.add( cons );
	cplex.extract(model);
}

void SUCsubprob::setMasterSoln ( vector< vector<bool> > & gen_stat ) {
	this->gen_stat = &gen_stat;
	
	int c=0;
	
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);

		for (int t=0; t<numPeriods; t++, c++) {
            if (genPtr->isMustUse) {
                cons[c].setBounds( gen_stat[g][t] * genPtr->maxCapacity, gen_stat[g][t] * genPtr->maxCapacity);
            } else {
                cons[c].setUB( gen_stat[g][t] * genPtr->maxCapacity );
            }
		}
	}

	// production amounts
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++, c++) {
			cons[c].setLB( gen_stat[g][t] * minGenerationReq[g] );
		}
	}
    
    // rest of the constraints are not a function of x
}

bool SUCsubprob::solve() {

	cutCoefs.reset();
	recourse_obj_val = 0;
	
	bool status = true;
	
	// process second-stage scenario subproblems
	for (int s=0; s<numScen; s++)
	{
		setup_subproblem (s);			// prepare the subproblem
		
		status = cplex.solve();			// solve the subproblem
		if (!status) {
			if ( cplex.getCplexStatus() == IloCplex::Infeasible ) {
				get_feasibility_cut_coefs(s);
				recourse_obj_val = INFINITY;
			}
			else {
				cout << "!! Error: Cplex returned the status " << cplex.getCplexStatus() << endl;
			}
			
			return false;
		}
		
		recourse_obj_val += sceProb[s] * cplex.getObjValue();	// get the obj value
		
		update_optimality_cut_coefs(s);	// update benders cut coefs
	}
	
	return true;
}

void SUCsubprob::setup_subproblem(int &s)
{
	// iterate over capacity constraints
	int c=0;
	
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		auto it = scenSet->mapVarNamesToIndex.find(genPtr->name);
		if ( it != scenSet->mapVarNamesToIndex.end() ) {		// if random supply generator
			
			for (int t=0; t<numPeriods; t++, c++) {				// Note: c is iterated in the secondary-loops
				if ((*gen_stat)[g][t]) {
					if ( genPtr->isMustUse ) {
						cons[c].setBounds( scenSet->vals[s][t*numBaseTimePerPeriod][it->second], scenSet->vals[s][t*numBaseTimePerPeriod][it->second]);
					}
					else {
						cons[c].setUB( scenSet->vals[s][t*numBaseTimePerPeriod][it->second] );
					}
				}
			}
		} else {
			c += numPeriods;
		}
	}
}

void SUCsubprob::update_optimality_cut_coefs(int &s)
{
	// get dual multipliers
	cplex.getDuals(duals, cons);
	
	// iterate over the constraints to write the Benders' cut
	int c=0;
	
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		auto it = scenSet->mapVarNamesToIndex.find(genPtr->name);
		if ( it != scenSet->mapVarNamesToIndex.end() ) {		// if random supply generator
			for (int t=0; t<numPeriods; t++, c++) {				// Note: c is iterated in the secondary-loops
				cutCoefs.pi_T[g][t] += sceProb[s] * duals[c] * ( scenSet->vals[s][t*numBaseTimePerPeriod][it->second] );
			}
		}
		else {
			for (int t=0; t<numPeriods; t++, c++) {				// Note: c is iterated in the secondary-loops
				cutCoefs.pi_T[g][t] += sceProb[s] * duals[c] * ( genPtr->maxCapacity );
			}
		}
	}
	
	// production amounts
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++, c++) {
			cutCoefs.pi_T[g][t] += sceProb[s] * duals[c] * minGenerationReq[g];
		}
	}
	
	// ramp up constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += sceProb[s] * duals[c] * genPtr->rampUpLim * periodLength;
		}
	}
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += sceProb[s] * duals[c] * genPtr->rampDownLim * periodLength;
		}
	}
	
	// model-specific constraints
	if (modelType == System) {
		// aggregated demand constraints
		for (int t=0; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += sceProb[s] * duals[c] * sysLoad[t];
		}
	}
	else {
		// Flow upper bounds
		for (int l=0; l<numLine; l++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += sceProb[s] * duals[c] * inst->powSys->lines[l].maxFlowLim;
			}
		}

		// Flow lower bounds
		for (int l=0; l<numLine; l++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += sceProb[s] * duals[c] * inst->powSys->lines[l].minFlowLim;
			}
		}
		
		// Phase-angle upper bounds
		for (int b=0; b<numBus; b++) {
			Bus *busPtr = &(inst->powSys->buses[b]);
			
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += sceProb[s] * duals[c] * ( busPtr->maxPhaseAngle - busPtr->minPhaseAngle );
			}
		}

		// DC-approximation to AC flow
		// - skipped, as all necessary coefs are 0
		// - important: these constraints are not added to cons array

		// Flow balance
		for (int b=0; b<numBus; b++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += sceProb[s] * duals[c] * busLoad[b][t];
			}
		}
	}
}

double SUCsubprob::getRecourseObjValue() {
	return recourse_obj_val;
}

void SUCsubprob::get_feasibility_cut_coefs(int &s)
{
	// reset earlier (optimality-)cut coefficient entries
	cutCoefs.reset();
	
	// get the extreme ray (based on the algorithm used to solve the subproblem)
	if ( cplex.getParam(IloCplex::RootAlg) == IloCplex::Dual ) {
	
		// Note: dualFarkas function mixes up the order of constraints
		//       We initialize an empty farkas constraint set, an empty dual-multiplier
		//		 vector, and a map to retrieve the correct values.
		IloRangeArray	farkasCons (env);
		
		// get the ray
		cplex.dualFarkas(farkasCons, duals);		// (duals is the ray...)
		
		// initialize the map
		farkasMap.clear();
		for (int c=0; c<farkasCons.getSize(); c++) {
			farkasMap.insert( pair<IloInt, double> (farkasCons[c].getId(), duals[c]) );
		}
	}
	else {
		cout << "!!Error: Cannot retrieve the extreme ray" << endl;
		exit(1);
	}
	
	// iterate over the constraints to write the Benders' feasibility cut
	int c=0;
	
	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		auto it = scenSet->mapVarNamesToIndex.find(genPtr->name);
		if ( it != scenSet->mapVarNamesToIndex.end() ) {		// if random supply generator
			for (int t=0; t<numPeriods; t++, c++) {				// Note: c is iterated in the secondary-loops
				cutCoefs.pi_T[g][t] += farkasMap[ cons[c].getId() ] * ( scenSet->vals[s][t*numBaseTimePerPeriod][it->second] );
			}
		}
		else {
			for (int t=0; t<numPeriods; t++, c++) {				// Note: c is iterated in the secondary-loops
				cutCoefs.pi_T[g][t] += farkasMap[ cons[c].getId() ] * ( genPtr->maxCapacity );
			}
		}
	}
	
	// production amounts
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++, c++) {
			cutCoefs.pi_T[g][t] += farkasMap[ cons[c].getId() ] * minGenerationReq[g];
		}
	}
	
	// ramp up constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * genPtr->rampUpLim * periodLength;
		}
	}
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=1; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * genPtr->rampDownLim * periodLength;
		}
	}
	
	// model-specific constraints
	if (modelType == System) {
		// aggregated demand constraints
		for (int t=0; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * sysLoad[t];
		}
	}
	else {
		// Flow upper bounds
		for (int l=0; l<numLine; l++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * inst->powSys->lines[l].maxFlowLim;
			}
		}
		
		// Flow lower bounds
		for (int l=0; l<numLine; l++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * inst->powSys->lines[l].minFlowLim;
			}
		}
		
		// Phase-angle upper bounds
		for (int b=0; b<numBus; b++) {
			Bus *busPtr = &(inst->powSys->buses[b]);
			
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * ( busPtr->maxPhaseAngle - busPtr->minPhaseAngle );
			}
		}
		
		// DC-approximation to AC flow
		// - skipped, as all necessary coefs are 0
		// - important: these constraints are not added to cons array
		
		// Flow balance
		for (int b=0; b<numBus; b++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * busLoad[b][t];
			}
		}
	}
}

/*
double SUCsubprob::computeLowerBound ()
{
    double bound = 0;

    IloEnv          env;
    IloModel        model   (env);
    IloRangeArray   cons    (env);
    IloCplex        cplex   (env);
    IloExpr         obj_func(env);

    /** Formulate the bounding problem **/

    /** Decision Variables **
    IloArray< IloNumVarArray > x (env, inst->G);    // generator states
    IloArray< IloNumVarArray > s (env, inst->G);    // start up
    IloArray< IloNumVarArray > z (env, inst->G);    // shut down
    IloArray< IloNumVarArray > p (env, inst->G);	// production amounts
    for (int g=0; g<inst->G; g++) {
        x[g] = IloNumVarArray(env, inst->T, 0, 1, ILOFLOAT);
        s[g] = IloNumVarArray(env, inst->T, 0, 1, ILOFLOAT);
        z[g] = IloNumVarArray(env, inst->T, 0, 1, ILOFLOAT);
        p[g] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
    }
    
    /** Constraints **
    // state constraints
    for (int g=0; g<inst->G; g++) {
        
        // t=0: generators are assumed to be turned on
        model.add( x[g][0] - 1 == s[g][0] - z[g][0] );
        
        // t>0
        for (int t=1; t<inst->T; t++) {
            model.add( x[g][t] - x[g][t-1] == s[g][t] - z[g][t] );
        }
    }
    
    // minimum uptime/downtime constraints
    for (int g=0; g<inst->G; g++)
    {
        // turn on inequalities
        for (int t=1; t<=inst->T; t++)
        {
            IloExpr lhs (env);
            for (int i = t-inst->min_u_time[g]+1; i<=t; i++)
            {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += 0;	// otherwise, generator is assumed to be operational (but turned on way earlier in the past)
            }
            model.add( lhs <= x[g][t-1] );
            lhs.end();
        }
    }
    
    for (int g=0; g<inst->G; g++)
    {
        // turn off inequalities
        for (int t=1; t<=inst->T; t++)
        {
            IloExpr lhs (env);
            for (int i = t-inst->min_d_time[g]+1; i<=t; i++)  {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += 0;	// otherwise, generator is assumed to be operational (but turned on way earlier in the past)
            }
            
            if (t-inst->min_d_time[g]-1 >= 0)	model.add( lhs <= 1 - x[g][t-inst->min_d_time[g]-1] );
            else								model.add( lhs <= 1 - 1 );	// assumed to be remaining on for a long, long, while
            lhs.end();
        }
    }

    // capacity constraints
    for (int g=0; g<inst->G; g++) {
        for (int t=0; t<inst->T; t++) {
            IloRange con;
            if (inst->must_run[g])  { con = IloRange (env, 0, p[g][t] - x[g][t] * inst->capacity[g], 0); }
            else                    { con = IloRange (env, -IloInfinity, p[g][t] - x[g][t] * inst->capacity[g], 0); }
            cons.add( con );
        }
    }

    // minimum production constraints
    for (int g=0; g<inst->G; g++) {
        for (int t=0; t<inst->T; t++) {
            model.add( p[g][t] >= inst->min_prod_lim[g] * x[g][t] );
        }
    }
    
    // ramping constraints
    for (int g=0; g<inst->G; g++) {
        for (int t=1; t<inst->T; t++) {
            model.add( p[g][t] - p[g][t-1] <= inst->ramp_u_lim[g] );
            model.add( p[g][t-1] - p[g][t] <= inst->ramp_d_lim[g] );
        }
    }
    
    // model-specific constraints, and the objective functions
    if (model_id == 0)
    {
        IloNumVarArray L (env, inst->T, 0, IloInfinity, ILOFLOAT);	// load shedding
        
        // aggregated demand constraints
        for (int t=0; t<inst->T; t++) {
            IloExpr expr (env);
            for (int g=0; g<inst->G; g++)	expr += p[g][t];
            expr += L[t];
            model.add( expr >= inst->aggDemand[t] );
        }
        
        // objective function
        for (int g=0; g<inst->G; g++) {
            for (int t=0; t<inst->T; t++) {
                obj_func += inst->var_cost[g] * p[g][t];
            }
        }
        
        for (int t=0; t<inst->T; t++) {
            obj_func += inst->load_shedding_penalty * L[t];
        }
    }
    else if (model_id == 1)
    {
        // node-based formulation's decision vars
        IloArray< IloNumVarArray > L (env, inst->B);        // load-shedding
        IloArray< IloNumVarArray > T (env, inst->B);		// phase angles
        for (int b=0; b<inst->B; b++) {
            L[b] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
            T[b] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
        }
        
        IloArray< IloNumVarArray > F (env, inst->L);		// flows
        for (int l=0; l<inst->L; l++) {
            F[l] = IloNumVarArray(env, inst->T, -IloInfinity, IloInfinity, ILOFLOAT);
        }
        
        // Flow upper and lower bounds
        for (int l=0; l<inst->L; l++) {
            for (int t=0; t<inst->T; t++) {
                model.add( inst->min_arc_flow[l] <= F[l][t] <= inst->max_arc_flow[l] );
            }
        }
        
        // Phase-angle upper bounds
        for (int b=0; b<inst->B; b++) {
            for (int t=0; t<inst->T; t++) {
                model.add( T[b][t] <= 2*inst->pi );
            }
        }
        
        // DC-approximation to AC flow
        for (int l=0; l<inst->L; l++) {
            for (int t=0; t<inst->T; t++) {
                int orig = inst->arcs[l].first;
                int dest = inst->arcs[l].second;
                
                model.add( F[l][t] == inst->susceptance[l] * (T[orig][t] - T[dest][t]) );
            }
        }
        
        // Flow balance
        for (int b=0; b<inst->B; b++) {
            for (int t=0; t<inst->T; t++) {
                
                IloExpr expr (env);
                
                // production
                for (int g=0; g<inst->G; g++) {
                    if (inst->generator_loc[g] == b) expr += p[g][t];
                }
                
                // incoming arcs
                for (int l=0; l<inst->L; l++) {
                    if (inst->arcs[l].second == b) expr += F[l][t];
                }
                
                // outgoing arcs
                for (int l=0; l<inst->L; l++) {
                    if (inst->arcs[l].first == b) expr -= F[l][t];
                }
                
                // load shedding
                expr += L[b][t];
                
                // constraint
                model.add( expr == inst->demand[b][t] );

                // free up memory
                expr.end();
            }
        }
        
        // objective function
        for (int g=0; g<inst->G; g++) {
            for (int t=0; t<inst->T; t++) {
                obj_func += inst->var_cost[g] * p[g][t];
            }
        }
        
        for (int b=0; b<inst->B; b++) {
            for (int t=0; t<inst->T; t++) {
                obj_func += inst->load_shedding_penalty * L[b][t];
            }
        }
    }
    else {
        cout << "!! Error: Unknown model id" << endl;
    }

    model.add( IloMinimize(env, obj_func) );
    model.add( cons );
    
    cplex.extract(model);
	/* TODO: Commented by HG after encountering "Invalid arguments candidates are: void setOut(? &)
    cplex.setOut(env.getNullStream()); */
    
    /** Change model components according to scenario, and solve the subproblems **
    for (int s=0; s<inst->S; s++) {
        
        // setup subproblem s
        for (int g=0, c=0; g<inst->G; g++) {    // only the capacity constraints contain randomness
            
            map<int, vector<vector<double> > >::iterator it = inst->rndSupply.find(g);
            
            if ( it != inst->rndSupply.end() ) {
                for (int t=0; t<inst->T; t++, c++) {	// Note: c is iterated in the secondary-loop
                    cons[c].setLinearCoef(x[g][t], -inst->rndSupply[g][s][t]);
                }
            }
            else {
                c += inst->T;
            }
        }    

        // solve the subproblem
        bool status = cplex.solve();
        if (!status) {
            cout << "Subproblem " << s << " returned the status " << cplex.getCplexStatus() << endl;
            exit(1);
        }
        
        bound += inst->sceProb[s] * cplex.getObjValue();
    }
    env.end();
    
    cout << "Computed recourse lower bound = " << bound << endl;
    
    return bound;
}
	 */
