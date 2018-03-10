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

	cons		= IloRangeArray (env);
	duals		= IloNumArray (env);

	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setParam(IloCplex::Threads, LShapedSubprobCPXThreads);
	cplex.setParam(IloCplex::PreInd, 0);
	//cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);	// Important
	
	solve_t = 0;
	cut_t	= 0;
	setup_t = 0;
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
	minGenerationReq.resize(numGen);					// minimum production requirements
	maxCapacity.resize(numGen);							// maximum generator capacities
	sysLoad.resize(numPeriods);							// aggregated system load
	resize_matrix(busLoad, numBus, numPeriods);			// individual bus loads
	
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

double SUCsubprob::getRandomCoef (int &s, int &t, int &loc) {
	if (t == 0 && probType == ShortTerm) {
		return inst->observations["RT"].vals[rep][(beginMin/periodLength)+(t*numBaseTimePerPeriod)][loc];
	}
	else {
		if (s == -1) {
			return (*expCapacity)[loc][t];
		} else {
			return inst->simulations.vals[s][(beginMin/periodLength)+(t*numBaseTimePerPeriod)][loc];
		}
	}
}

void SUCsubprob::formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep, IloArray<IloNumArray> &masterSoln, vector<vector<double>> &expCapacity)
{
	// get a handle on the instance and model type
	this->inst		= &inst;
	this->beginMin	= beginMin;
	this->probType	= probType;
	this->modelType = modelType;
	this->rep		= rep;
	this->genState	= &masterSoln;
	this->expCapacity = &expCapacity;
	
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
	x = IloArray< IloNumVarArray > (env, numGen);	// generator states
	for (int g=0; g<numGen; g++) {
		p[g] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
		x[g] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
		
		sprintf(buffer, "p_%d", g);
		p[g].setNames(buffer);
		
		sprintf(buffer, "x_%d", g);
		x[g].setNames(buffer);
	}

	/** Constraints **/

	// state constraints
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++) {
			cons.add( IloRange (env, 1, x[g][t], 1) );
		}
	}
	
	// capacity constraints
	// Variant 1: Leads to random coefficient matrix
//	for (int g=0; g<numGen; g++) {
//		Generator *genPtr = &(inst->powSys->generators[g]);
//
//		for (int t=0; t<numPeriods; t++) {
//			if (genPtr->isMustUse) {
//				cons.add( IloRange(env, 0, p[g][t] - x[g][t] * genPtr->maxCapacity, 0) );
//			} else {
//				cons.add( IloRange(env, -IloInfinity, p[g][t] - x[g][t] * genPtr->maxCapacity, 0) );
//			}
//		}
//	}
	// Variant 2: Leads to random rhs vector
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		for (int t=0; t<numPeriods; t++) {
			if (genPtr->isMustUse) {
				cons.add( IloRange(env, genPtr->maxCapacity, p[g][t], genPtr->maxCapacity) );
			} else {
				cons.add( IloRange(env, -IloInfinity, p[g][t], genPtr->maxCapacity) );
			}
		}
	}

	
	// minimum production constraints
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++) {
			model.add( p[g][t] >= x[g][t] * minGenerationReq[g] );
		}
	}
	
	// ramp up constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		/* t = 0 */
		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {	// prev generation is available
			IloRange con (env, -IloInfinity, p[g][t] - getGenProd(g, t-1), genPtr->rampUpLim * periodLength);
			cons.add(con);
		}
		
		/* t > 0 */
		for (t=1; t<numPeriods; t++) {
			IloRange con (env, -IloInfinity, p[g][t] - p[g][t-1], genPtr->rampUpLim * periodLength);
			cons.add(con);
		}
	}	
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		// input-inconsistency check
		int shutDownPeriod = checkShutDownRampDownInconsistency(g);

		double rampDownRate;
		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			IloRange con (env, -IloInfinity, getGenProd(g, t-1) - p[g][t], rampDownRate);
			cons.add(con);
		}
		for (t=1; t<numPeriods; t++) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			IloRange con (env, -IloInfinity, p[g][t-1] - p[g][t], rampDownRate);
			cons.add(con);
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
	
	// No load-shedding in buses with 0 load
	for (int b=0; b<numBus; b++) {
		for (int t=0; t<numPeriods; t++) {
			if (busLoad[b][t] < EPSzero) L[b][t].setUB(0);
		}
	}
	
	// No over-generation in buses with no generators
	for (int b=0; b<numBus; b++) {
		if (inst->powSys->buses[b].connectedGenerators.size() == 0) {
			for (int t=0; t<numPeriods; t++) {
				O[b][t].setUB(0);
			}
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
			for (int g=0; g < (int) busPtr->connectedGenerators.size(); g++) {
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
			obj -= genPtr->variableCost * periodLength/60.0 * minGenerationReq[g] * x[g][t];
		}
	}
	
	for (int b=0; b<numBus; b++) {
		for (int t=0; t<numPeriods; t++) {
			obj += loadShedPenaltyCoef * L[b][t];
			obj += overGenPenaltyCoef * O[b][t];
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
		for (int g=0; g<numGen; g++) expr += p[g][t];
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
		obj += overGenPenaltyCoef * O[t];
	}
	
	// prepare the model and the solver
	model.add( IloMinimize(env, obj) );
	model.add( cons );
	cplex.extract(model);
}

void SUCsubprob::setMasterSoln () {
	int c=0;

	IloRangeArray stateCons (env);
	IloNumArray stateVals (env);
	
	// state constraints
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++) {
			stateCons.add(cons[c]);
			stateVals.add( (double)(*genState)[g][t] );
			c++;
		}
	}
	stateCons.setBounds(stateVals, stateVals);
	
	stateCons.end();
	stateVals.end();

	// rest of the constraints are not a function of x
}

bool SUCsubprob::solve(int mappedScen, BendersCutCoefs &cutCoefs, double &objValue, vector<double> &initGen) {
	
	cutCoefs.reset();
	
	setup_subproblem (mappedScen);			// prepare the subproblem
	bool status = cplex.solve();	// solve the subproblem
	
	// feasibility
	if (status) {
		compute_optimality_cut_coefs(cutCoefs, mappedScen);
		objValue = cplex.getObjValue();
		
		if (probType==DayAhead || (probType==ShortTerm && beginMin==0))	getInitGen(initGen);
	}
	
	// infeasibility
	if (cplex.getCplexStatus() == IloCplex::Infeasible ) {
		compute_feasibility_cut_coefs(cutCoefs, mappedScen);
		objValue = INFINITY;
	}
	
	return status;
}

void SUCsubprob::getInitGen(vector<double> &initGen) {
	for (int g=0; g<numGen; g++) {
		initGen[g] = cplex.getValue( p[g][0] );
	}
}

void SUCsubprob::setup_subproblem(int &s) {
	double supply;
	int c=numGen*numPeriods;

	// capacity constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		auto it = inst->simulations.mapVarNamesToIndex.find(genPtr->name);
		
		if ( it != inst->simulations.mapVarNamesToIndex.end() ) {	// if random supply generator
			for (int t=0; t<numPeriods; t++, c++) {					// Note: c is iterated in the secondary-loops
				// Variant 1: Modifying random coefficient matrices
//				if (s==-1)	{ supply = min(getRandomCoef(s, t, g), genPtr->maxCapacity); }
//				else 		{ supply = min(getRandomCoef(s, t, it->second), genPtr->maxCapacity); }
//
//				cons[c].setLinearCoef(x[g][t], -supply);

				// Variant 2: Modifying random rhs vector
				if ( fabs((double)(*genState)[g][t]) < 0.5 ) {
					supply = 0;
				} else {
					if (s==-1)	{ supply = min(getRandomCoef(s, t, g), genPtr->maxCapacity); }
					else 		{ supply = min(getRandomCoef(s, t, it->second), genPtr->maxCapacity); }
				}
				
				if (genPtr->isMustUse) {
					cons[c].setBounds(supply, supply);	// ax = b
				} else {
					cons[c].setUB(supply);	// ax <= b
				}
			}
		} else {
//			// Variant 1:
//			c += numPeriods;
			// Variant 2:
			for (int t=0; t<numPeriods; t++, c++) {					// Note: c is iterated in the secondary-loops
				if ( fabs((double)(*genState)[g][t]) < 0.5 ) {
					supply = 0;
				} else {
					supply = genPtr->maxCapacity;
				}
				if (genPtr->isMustUse) {
					cons[c].setBounds(supply, supply);	// ax = b
				} else {
					cons[c].setUB(supply);	// ax <= b
				}
			}
		}
	}
}

void SUCsubprob::compute_optimality_cut_coefs(BendersCutCoefs &cutCoefs, int &s)
{
	// get dual multipliers
	cplex.getDuals(duals, cons);
	
	// iterate over the constraints to write the Benders' cut
	int c=0;
	
	// state constraints
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++) {
			cutCoefs.pi_T[g][t] += duals[c];
			c++;
		}
	}
	
	// capacity
	// Variant 1: Randomness is in the coefficient matrix
//	// skipped, as all necessary coefs are 0
//	c += numGen * numPeriods;
	
	// Variant 2: Randomness is in the rhs vector
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		double supply;
		
		auto it = inst->simulations.mapVarNamesToIndex.find(genPtr->name);
		
		if ( it != inst->simulations.mapVarNamesToIndex.end() ) {	// if random supply generator
			for (int t=0; t<numPeriods; t++, c++) {					// Note: c is iterated in the secondary-loops
				if (s==-1)	{ supply = min(getRandomCoef(s, t, g), genPtr->maxCapacity); }
				else 		{ supply = min(getRandomCoef(s, t, it->second), genPtr->maxCapacity); }
				cutCoefs.pi_T[g][t] += duals[c] * supply;
			}
		} else {
			for (int t=0; t<numPeriods; t++, c++) {					// Note: c is iterated in the secondary-loops
				supply = genPtr->maxCapacity;
				cutCoefs.pi_T[g][t] += duals[c] * supply;
			}
		}
	}
	
	// minimum generation requirements
	// skipped, as all necessary coefs are 0
	
	// ramp up constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);

		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {	// prev generation is available
			cutCoefs.pi_b += duals[c] * (genPtr->rampUpLim*periodLength + getGenProd(g, t-1));
			c++;
		}
		for (t=1; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += duals[c] * genPtr->rampUpLim * periodLength;
		}
	}
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		// input-inconsistency check
		int shutDownPeriod = checkShutDownRampDownInconsistency(g);
		
		double rampDownRate;
		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
            
			cutCoefs.pi_b += duals[c] * (rampDownRate - getGenProd(g,t-1));
			c++;
		}
		for (t=1; t<numPeriods; t++, c++) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;

			cutCoefs.pi_b += duals[c] * rampDownRate;
		}
	}
	
	// model-specific constraints
	if (modelType == System) {
		// aggregated demand constraints
		for (int t=0; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += duals[c] * sysLoad[t];
		}
	}
	else {
		// Flow upper bounds
		for (int l=0; l<numLine; l++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += duals[c] * inst->powSys->lines[l].maxFlowLim;
			}
		}

		// Flow lower bounds
		for (int l=0; l<numLine; l++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += duals[c] * inst->powSys->lines[l].minFlowLim;
			}
		}
		
		// Phase-angle upper bounds
		for (int b=0; b<numBus; b++) {
			Bus *busPtr = &(inst->powSys->buses[b]);
			
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += duals[c] * ( busPtr->maxPhaseAngle - busPtr->minPhaseAngle );
			}
		}
		
		// No load-shedding in buses with 0 load
		// - skipped, as all necessary coefs are 0

		// No over-generation in buses with no generators
		// - skipped, as all necessary coefs are 0

		// DC-approximation to AC flow
		// - skipped, as all necessary coefs are 0
		// - important: these constraints are not added to cons array

		// Flow balance
		for (int b=0; b<numBus; b++) {
			for (int t=0; t<numPeriods; t++, c++) {
				cutCoefs.pi_b += duals[c] * busLoad[b][t];
			}
		}
	}
}

void SUCsubprob::compute_feasibility_cut_coefs(BendersCutCoefs &cutCoefs, int &s)
{
	// reset earlier (optimality-)cut coefficient entries
	cutCoefs.reset();
	
	// resolve the problem with presolve turned off (oth. basis may not be available)
	cplex.setParam(IloCplex::PreInd, 0);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
	cplex.solve();
	cplex.setParam(IloCplex::RootAlg, IloCplex::Auto);
	cplex.setParam(IloCplex::PreInd, 1);
	
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
	
	// state constraints
	for (int g=0; g<numGen; g++) {
		for (int t=0; t<numPeriods; t++) {
			cutCoefs.pi_T[g][t] += farkasMap[ cons[c].getId() ];
			c++;
		}
	}
	
	// capacity constraints
	// Variant 1: Randomness is in the coefficient matrix
//	// skipped, as all necessary coefs are 0
//	c += numGen*numPeriods;
	// Variant 2: Randomness is in the rhs vector
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		double supply;
		
		auto it = inst->simulations.mapVarNamesToIndex.find(genPtr->name);
		
		if ( it != inst->simulations.mapVarNamesToIndex.end() ) {	// if random supply generator
			for (int t=0; t<numPeriods; t++, c++) {					// Note: c is iterated in the secondary-loops
				if (s==-1)	{ supply = min(getRandomCoef(s, t, g), genPtr->maxCapacity); }
				else 		{ supply = min(getRandomCoef(s, t, it->second), genPtr->maxCapacity); }
				cutCoefs.pi_T[g][t] += farkasMap[ cons[c].getId() ] * supply;
			}
		} else {
			for (int t=0; t<numPeriods; t++, c++) {					// Note: c is iterated in the secondary-loops
				supply = genPtr->maxCapacity;
				cutCoefs.pi_T[g][t] += farkasMap[ cons[c].getId() ] * supply;
			}
		}
	}
	
	// ramp up constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {	// prev generation is available
			cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * (genPtr->rampUpLim*periodLength + getGenProd(g, t-1));
			c++;
		}
		for (t=1; t<numPeriods; t++, c++) {
			cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * genPtr->rampUpLim * periodLength;
		}
	}
	
	// ramp down constraints
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		// input-inconsistency check
		int shutDownPeriod = checkShutDownRampDownInconsistency(g);
		
		double rampDownRate;
		int t=0;
		if (getGenProd(g, t-1) > -INFINITY) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			cutCoefs.pi_b += farkasMap[ cons[c].getId() ] * (rampDownRate - getGenProd(g,t-1));
			c++;
		}
		for (t=1; t<numPeriods; t++, c++) {
			rampDownRate = (shutDownPeriod != t) ? genPtr->rampDownLim*periodLength : genPtr->maxCapacity;
			
			cutCoefs.pi_b += farkasMap[ cons[c].getId() ]  * rampDownRate;
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


/****************************************************************************
 * getEDGenProd
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded generation of the generator by the ED model.
 ****************************************************************************/
double SUCsubprob::getEDGenProd(int genId, int period) {
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
 * getGenState
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded status of the generator.
 ****************************************************************************/
bool SUCsubprob::getGenState(int genId, int period) {
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
 * getGenProd
 * - If period != the beginning of the planning horizon, returns getEDGenProd
 * - If period is the beginning of the planning horizon:
 * 		- Use generator histories from the previous run:
 *		- Don't use generator histories
 *			- DA-UC: returns -infinity due to lack of info
 *			- ST-UC: returns 1st-period generation amount from the DA solve,
 * 			if the gen is DA, -infinity if it is ST.
 ****************************************************************************/
double SUCsubprob::getGenProd(int g, int t) {
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
int SUCsubprob::checkShutDownRampDownInconsistency (int g) {
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
double SUCsubprob::getUCGenProd(int genId, int period) {
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