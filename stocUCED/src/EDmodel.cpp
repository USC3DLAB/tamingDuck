/*
 * EDmodel.cpp
 *
 *  Created on: Dec 12, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "EDmodel.hpp"

extern runType runParam;

EDmodel::EDmodel() {
	model = IloModel(env);
	cplex = IloCplex(env);
}

EDmodel::~EDmodel() {
	env.end();
}

void EDmodel::preprocess(instance &inst, int t0) {

	/* Parameters to be updated */
	resize_matrix(minGenerationReq, inst.powSys->numGen, runParam.ED_numPeriods);
	resize_matrix(maxGenerationReq, inst.powSys->numGen, runParam.ED_numPeriods);
	resize_matrix(busLoad, inst.powSys->numBus, runParam.ED_numSolves);

	for (int g = 0; g < inst.powSys->numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);

		// Make sure that the generators meet the minGenerationReq & Assumption 1
		double temp = genPtr->minGenerationReq;
		if (temp > min(genPtr->rampUpLim * runParam.ED_resolution, genPtr->rampDownLim * runParam.ED_resolution)) {
			temp = min(genPtr->rampUpLim * runParam.ED_resolution, genPtr->rampDownLim * runParam.ED_resolution);
		}
		fill(minGenerationReq[g].begin(), minGenerationReq[g].begin()+runParam.ED_numPeriods, temp);
		fill(maxGenerationReq[g].begin(), maxGenerationReq[g].begin()+runParam.ED_numPeriods, genPtr->maxCapacity);

		/* Make sure that the generator minimum is set to 0. */
		for ( int t = t0; t < t0 + runParam.ED_numPeriods; t++ ) {
			minGenerationReq[g][t-t0] *= (double) (inst.solution.x[g][t] && inst.solution.s[g][t]);
			maxGenerationReq[g][t-t0] *= (double) (inst.solution.x[g][t] && inst.solution.s[g][t]);
		}
	}

	/* Compute the load at a bus based on the regional load */
	for (int b = 0; b < inst.powSys->numBus; b++) {
		Bus *busPtr = &(inst.powSys->buses[b]);

		// TODO: Fix the mapping issue
		// int r = inst.detObserv[2].mapVarNamesToIndex[numToStr(busPtr->regionId)];
		for (int t = 0; t < runParam.numPeriods; t++) {
			busLoad[b][t] = 10;//inst.detObserv[2].vals[0][t0+t][0]*busPtr->loadPercentage;
		}
	}


	return;
}

void EDmodel::formulate(instance &inst, int beginTime) {
	char elemName[NAMESIZE];

	/**** Decision variables *****/
	/* Generation and overgeneration */
	gen = IloArray< IloNumVarArray > (env, inst.powSys->numGen);
	overGen = IloArray< IloNumVarArray > (env, inst.powSys->numGen);
	for ( int g = 0; g < inst.powSys->numGen; g++ ) {

		sprintf(elemName, "gen[%d]", g);
		gen[g] = IloNumVarArray(env, runParam.ED_numPeriods, 0, IloInfinity, ILOFLOAT);
		gen[g].setNames(elemName); model.add(gen[g]);

		sprintf(elemName, "overGen[%d]", g);
		overGen[g] = IloNumVarArray(env, runParam.ED_numPeriods, 0, IloInfinity, ILOFLOAT);
		overGen[g].setNames(elemName); model.add(overGen[g]);
	}

	/* Demand met, demand shed and bus angle */
	demMet = IloArray< IloNumVarArray > (env, inst.powSys->numBus);
	demShed = IloArray <IloNumVarArray> (env, inst.powSys->numBus);
	theta = IloArray< IloNumVarArray > (env, inst.powSys->numBus);
	for ( int b = 0; b < inst.powSys->numBus; b++ ) {
		Bus *bptr = &(inst.powSys->buses[b]);

		sprintf(elemName, "demMet[%d]", b);
		demMet[b] = IloNumVarArray(env, runParam.ED_numPeriods, 0, IloInfinity, ILOFLOAT);
		demMet[b].setNames(elemName); model.add(demMet[b]);

		sprintf(elemName, "demShed[%d]", b);
		demShed[b] = IloNumVarArray(env, runParam.ED_numPeriods, 0, IloInfinity, ILOFLOAT);
		demShed[b].setNames(elemName); model.add(demShed[b]);

		sprintf(elemName, "theta[%d]", b);
		theta[b] = IloNumVarArray(env, runParam.ED_numPeriods, bptr->minPhaseAngle, bptr->maxPhaseAngle, ILOFLOAT);
		theta[b].setNames(elemName); model.add(theta[b]);
	}

	/* Power flows */
	flow = IloArray< IloNumVarArray > (env, inst.powSys->numLine);
	for ( int l = 0; l < inst.powSys->numLine; l++ ) {
		Line *lptr = &(inst.powSys->lines[l]);
		sprintf(elemName, "flow[%d,%d]", lptr->orig->id, lptr->dest->id);

		flow[l] = IloNumVarArray(env, runParam.ED_numPeriods, lptr->minFlowLim, lptr->maxFlowLim, ILOFLOAT);
		flow[l].setNames(elemName); model.add(flow[l]);
	}

	/***** Constraints *****/
	/* Flow balance equation */
	for (int t = 0; t < runParam.ED_numPeriods; t++ ) {
		for ( int b = 0; b < inst.powSys->numBus; b++ ) {
			IloExpr expr (env);
			sprintf(elemName, "flowBalance[%d][%d]", b, t);

			for ( int g = 0; g < inst.powSys->numGen; g++ )
				if ( inst.powSys->generators[g].connectedBus->id == b ) expr += (gen[g][t] - overGen[b][t]);
			for ( int l = 0; l < inst.powSys->numLine; l++ ) {
				if ( inst.powSys->lines[l].dest->id == b ) expr += flow[l][t];
				if ( inst.powSys->lines[l].orig->id == b ) expr -= flow[l][t];
			}
			expr -= demMet[b][t];

			IloConstraint c( expr == 0 ); c.setName(elemName); model.add(c);
			expr.end();
		}
	}

	/* Line power flow equation : DC approximation */
	for (int t = 0; t < runParam.ED_numPeriods; t++ ) {
		for ( int l = 0; l < inst.powSys->numLine; l++ ) {
			int orig = inst.powSys->lines[l].orig->id;
			int dest = inst.powSys->lines[l].dest->id;
			sprintf(elemName, "dcApprox[%d,%d][%d]", orig, dest, t);

			IloConstraint c(  flow[l][t] == inst.powSys->lines[l].susceptance*(theta[orig][t] - theta[dest][t]) ); c.setName(elemName); model.add(c);
		}
	}

	/* Generation ramping constraints */
	{
		int t = 0;
		for ( int g = 0; g < inst.powSys->numGen; g++ ) {
			if ( inst.solution.x[beginTime][g] == 1 ) {
				sprintf(elemName, "rampUp[%d][%d]", g, t);
				IloConstraint c1( gen[g][t] -  inst.solution.g[beginTime-1][g] <= inst.powSys->generators[g].rampUpLim); c1.setName(elemName); model.add(c1);
				sprintf(elemName, "rampDown[%d][%d]", g, t);
				IloConstraint c2( inst.powSys->generators[g].rampDownLim <= gen[g][t] - inst.solution.g[beginTime-1][g]); c2.setName(elemName); model.add(c2);
			}
		}
	}
	for ( int t = 1; t < runParam.ED_numPeriods; t++ ) {
		for ( int i = 0; i < inst.powSys->numGen; i++ ) {
			sprintf(elemName, "rampUp[%d][%d]", i, t);
			IloConstraint c1( gen[i][t] - gen[i][t-1] <= inst.powSys->generators[i].rampUpLim); c1.setName(elemName); model.add(c1);
			sprintf(elemName, "rampDown[%d][%d]", i, t);
			IloConstraint c2( inst.powSys->generators[i].rampDownLim <= gen[i][t] - gen[i][t-1]); c2.setName(elemName); model.add(c2);
		}
	}

	/* Generation capacity and availability limit */
	for ( int g = 0; g < inst.powSys->numGen; g++ ) {
		for (int t = 0; t < runParam.ED_numPeriods; t++ ) {
			/* Make sure that the generation capacity limits are met. The minimum and maximum for generators which are turned off are set to zero during pre-processing */
			gen[g][t].setBounds(minGenerationReq[g][t], maxGenerationReq[g][t]);

			/* Make sure that the over generation is less than the actual generation!! */
			sprintf(elemName, "genConsist[%d][%d]", g, t);
			IloConstraint c (overGen[g][t] - gen[g][t] <= 0); c.setName(elemName); model.add(c);
		}
	}

	/* Demand consistency */
	for ( int d = 0; d < inst.powSys->numBus; d++ ) {
		for ( int t = 0; t < runParam.ED_numPeriods; t++) {
			sprintf(elemName, "demConsist[%d][%d]", d, t);
			IloConstraint c(demMet[d][t] + demShed[d][t] == busLoad[d][t]); c.setName(elemName); model.add(c);
		}
	}

	/***** Objective function *****/
	IloExpr realTimeCost (env);
	IloObjective obj;
	for ( int t = 0; t < runParam.ED_numPeriods; t++) {
		/* Generation cost */
		for ( int g = 0; g < inst.powSys->numGen; g++ )
			realTimeCost += (inst.powSys->generators[g].variableCost*runParam.ED_resolution/60)*gen[g][t];

		/* Load shedding penalty */
		for ( int d = 0; d < inst.powSys->numBus; d++ )
			realTimeCost += loadShedPenaltyCoef*demShed[d][t];
	}
	obj = IloMinimize(env, realTimeCost);
	model.add(obj);
	realTimeCost.end();

	/* Model parameters */
	cplex.extract(model);

	cplex.exportModel("EDform.lp");

}//END formulate()

/****************************************************************************
 * solve
 * - Solves the model and records the solution into inst->solution.
 ****************************************************************************/
bool EDmodel::solve(instance &inst, int t0) {
	bool status;

	try {
		status = cplex.solve();

		// record the solution
		if (status) {
			for (int g = 0; g < inst.powSys->numGen; g++) {
				Generator *genPtr = &(inst.powSys->generators[g]);
				for (int t = 0; t < runParam.ED_numPeriods; t++) {
					inst.solution.g[g][t0+t] = cplex.getValue(gen[g][t]);
				}
			}
		}
	}
	catch (IloException &e) {
		cout << e << endl;
	}

	return status;
}
