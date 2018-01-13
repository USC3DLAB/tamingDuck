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

#include "misc.hpp"
#include "EDmodel.hpp"

extern runType runParam;

EDmodel::EDmodel(instance &inst, int t0, int rep) {

	model = IloModel(env);
	cplex = IloCplex(model);

	cplex.setOut(env.getNullStream());

	/* Pre-process to setup some parameters: generator status (ramping and capacity), load, and stochastic generation */
	numBus = numLoad = inst.powSys->numBus;
	numGen = inst.powSys->numGen;
	numLine = inst.powSys->numLine;

	numPeriods = runParam.ED_numPeriods;

	resize_matrix(busLoad, numBus, numPeriods);
	resize_matrix(genMin, numGen, numPeriods);
	resize_matrix(genMax, numGen, numPeriods);
	resize_matrix(genRampDown, numGen, numPeriods);
	resize_matrix(genRampUp, numGen, numPeriods);
	resize_matrix(genAvail, numGen, numPeriods);

	for ( int d = 0; d < numBus; d++ ) {
		Bus *busPtr = &(inst.powSys->buses[d]);
		int r = inst.detObserv[2].mapVarNamesToIndex[num2str(busPtr->regionId)];

		for ( int t = 0; t < runParam.ED_numPeriods; t++) {
			busLoad[d][t] = inst.detObserv[2].vals[rep][t0+t][r]*busPtr->loadPercentage;
		}
	}

	for ( int g = 0; g < numGen; g++ ) {
		Generator genPtr = inst.powSys->generators[g];
		for (int t = 0; t < numPeriods; t++ ) {
			/* If the model horizon exceeds the number of periods, then the last period commitment are repeated for the remainder of horizon */
			int idxT = min(t0+t, runParam.numPeriods);

			/* Make sure that the generation capacity limits are met. The minimum and maximum for generators which are turned off are set to zero during pre-processing */
			genMax[g][t] = inst.solution.x[g][idxT]*
					genPtr.maxCapacity;
			genMin[g][t] = inst.solution.x[g][idxT]*
					min(genPtr.minGenerationReq, min(genPtr.rampUpLim * runParam.ED_resolution, genPtr.rampDownLim * runParam.ED_resolution));

			/* Stochastic generation set to what is available */
			auto it = inst.stocObserv[2].mapVarNamesToIndex.find(genPtr.name);
			if ( it != inst.stocObserv[2].mapVarNamesToIndex.end() ) {
				genAvail[g][t] = min(genMax[g][t], inst.stocObserv[2].vals[rep][idxT][g]);
			}

			if ( idxT != 0 ) {
				/* It is not the first period (t0 == 0 and t == 0) of the entire run, so we have access to solutions */
				genRampUp[g][t] = inst.solution.x[g][idxT-1]*inst.powSys->generators[g].rampUpLim +
						inst.solution.x[g][idxT]*(1-inst.solution.x[g][idxT-1])*inst.powSys->generators[g].maxCapacity;
				genRampDown[g][t] = inst.solution.x[g][idxT]*inst.powSys->generators[g].rampDownLim +
						inst.solution.x[g][idxT-1]*(1-inst.solution.x[g][idxT])*inst.powSys->generators[g].maxCapacity;
			}
		}
	}

}

EDmodel::~EDmodel() {
	env.end();
}

void EDmodel::formulate(instance &inst, int t0) {
	char elemName[NAMESIZE];

	/* time decomposition and stoc elements */
	timeCols = vector <string> ();
	timeRows = vector <string> ();
	stocRows = vector <string> ();

	/**** Decision variables *****/
	genUsed = IloArray< IloNumVarArray > (env, numGen);
	overGen = IloArray< IloNumVarArray > (env, numGen);
	demMet = IloArray< IloNumVarArray > (env, numLoad);
	demShed = IloArray <IloNumVarArray> (env, numLoad);
	theta = IloArray< IloNumVarArray > (env, numLoad);
	flow = IloArray< IloNumVarArray > (env, numLine);

	for ( int t = 0; t < numPeriods; t++ ) {
		/* Generation and over-generation */
		for ( int g = 0; g < numGen; g++ ) {
			if ( t == 0 ) {
				genUsed[g] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
				overGen[g] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
			}

			sprintf(elemName, "genUsed[%d][%d]", g, t);
			genUsed[g][t].setName(elemName); model.add(genUsed[g][t]);

			if ( t == 0 && g == 0 )
				timeCols.push_back(elemName);
			else if ( t == 1 && g == 0 )
				timeCols.push_back(elemName);

			sprintf(elemName, "overGen[%d][%d]", g, t);
			overGen[g][t].setName(elemName); model.add(overGen[g][t]);
		}

		/* Demand met, demand shed and bus angle */
		for ( int b = 0; b < numBus; b++ ) {
			Bus *bptr = &(inst.powSys->buses[b]);

			if ( t == 0 ) {
				demMet[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
				demShed[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
				theta[b] = IloNumVarArray(env, numPeriods, bptr->minPhaseAngle, bptr->maxPhaseAngle, ILOFLOAT);
			}

			sprintf(elemName, "demMet[%d][%d]", b, t);
			demMet[b][t].setName(elemName); model.add(demMet[b][t]);

			sprintf(elemName, "demShed[%d][%d]", b, t);
			demShed[b][t].setName(elemName); model.add(demShed[b][t]);

			sprintf(elemName, "theta[%d][%d]", b, t);
			theta[b][t].setName(elemName); model.add(theta[b][t]);
		}

		/* Power flows */
		for ( int l = 0; l < numLine; l++ ) {
			Line *lptr = &(inst.powSys->lines[l]);

			if ( t == 0 ) {
				flow[l] = IloNumVarArray(env, numPeriods, lptr->minFlowLim, lptr->maxFlowLim, ILOFLOAT);
			}

			sprintf(elemName, "flow[%s][%d]", lptr->name.c_str(), t);
			flow[l][t].setName(elemName); model.add(flow[l][t]);
		}
	}

	/***** Constraints *****/
	for (int t = 0; t < numPeriods; t++ ) {
		/* Flow balance equation */
		for ( int b = 0; b < numBus; b++ ) {
			IloExpr expr (env);
			sprintf(elemName, "flowBalance[%d][%d]", b, t);

			for ( int g = 0; g < numGen; g++ )
				if ( inst.powSys->generators[g].connectedBus->id == b ) expr += genUsed[g][t];
			for ( int l = 0; l < numLine; l++ ) {
				if ( inst.powSys->lines[l].dest->id == b ) expr += flow[l][t];
				if ( inst.powSys->lines[l].orig->id == b ) expr -= flow[l][t];
			}
			expr -= demMet[b][t];

			IloConstraint c( expr == 0 ); c.setName(elemName); model.add(c);
			expr.end();

			if ( t == 0 && b == 0 )
				timeRows.push_back(elemName);
			else if ( t == 1 && b == 0 )
				timeRows.push_back(elemName);
		}

		/* Line power flow equation : DC approximation */
		for ( int l = 0; l < numLine; l++ ) {
			int orig = inst.powSys->lines[l].orig->id;
			int dest = inst.powSys->lines[l].dest->id;
			sprintf(elemName, "dcApprox[%s][%d]", inst.powSys->lines[l].name.c_str(), t);

			IloExpr expr (env);
			expr = flow[l][t] - inst.powSys->lines[l].susceptance*(theta[orig][t] - theta[dest][t]);

			IloConstraint c( expr == 0); c.setName(elemName); model.add(c);
		}

		/* Generation ramping constraints: applicable only if the generator continues to be ON, i.e., x[g][t] = 1. */
		if ( t == 0 ) {
			if ( t0 == 0 ) {
				/* TODO: The initial generation level is not considered */
				cout << "Warning: Initial generation ignored." << endl;
			}
			else {
				/* The first time period of the ED horizon */
				for ( int g = 0; g < numGen; g++ ) {
					/* ramp-up */
					sprintf(elemName, "rampUp[%d][%d]", g, t);
					IloConstraint c1( genUsed[g][t] + overGen[g][t] -  inst.solution.g_ED[g][t0-1] <= genRampUp[g][t]);
					c1.setName(elemName); model.add(c1);

					/* ramp-down */
					sprintf(elemName, "rampDown[%d][%d]", g, t);
					IloConstraint c2( inst.solution.g_ED[g][t0-1] - genUsed[g][t] - overGen[g][t] <= genRampDown[g][t]);
					c2.setName(elemName); model.add(c2);
				}
			}
		}
		else {
			/* The the remaining periods of the ED horizon */
			for ( int g = 0; g < numGen; g++ ) {
				/* ramp-up */
				sprintf(elemName, "rampUp[%d][%d]", g, t);
				IloConstraint c1( genUsed[g][t] + overGen[g][t] - genUsed[g][t-1] - overGen[g][t-1] <= genRampUp[g][t]);
				c1.setName(elemName); model.add(c1);

				/* ramp-down */
				sprintf(elemName, "rampDown[%d][%d]", g, t);
				IloConstraint c2( genUsed[g][t-1] + overGen[g][t-1] - genUsed[g][t] - overGen[g][t] <= genRampDown[g][t]);
				c2.setName(elemName); model.add(c2);
			}
		}

		/* TODO: Fix generation capacity for stochastic generators. Generation capacity and availability limit */
		for ( int g = 0; g < numGen; g++ ) {
			Generator genPtr = inst.powSys->generators[g];
			/* Make sure that the generation capacity limits are met.
			 * The minimum and maximum for generators which are turned off are set to zero during pre-processing */
			sprintf(elemName, "maxGen[%d][%d]", g, t);
			IloConstraint c1( genUsed[g][t] + overGen[g][t] <= genMax[g][t]);
			c1.setName(elemName); model.add(c1);

			sprintf(elemName, "minGen[%d][%d]", g, t);
			IloConstraint c2( genUsed[g][t] + overGen[g][t] >= genMin[g][t]);
			c2.setName(elemName); model.add(c2);

			/* Stochastic generation consistency */
			auto it = inst.stocObserv[2].mapVarNamesToIndex.find(genPtr.name);
			if ( it != inst.stocObserv[2].mapVarNamesToIndex.end() ) {
				sprintf(elemName, "stocAvail[%d][%d]", g, t);
				IloConstraint c (genUsed[g][t] + overGen[g][t] <= genAvail[g][t]); c.setName(elemName); model.add(c);

				if ( t != 0 )
					stocRows.push_back(elemName);
			}
		}

		/* Demand consistency */
		for ( int d = 0; d < numBus; d++ ) {
			sprintf(elemName, "demConsist[%d][%d]", d, t);
			IloConstraint c(demMet[d][t] + demShed[d][t] == busLoad[d][t]); c.setName(elemName); model.add(c);
		}
	}

	/***** Objective function *****/
	IloExpr realTimeCost (env);
	IloObjective obj;
	sprintf(elemName, "realTimeCost");
	for ( int t = 0; t < numPeriods; t++) {
		/* Generation cost */
		for ( int g = 0; g < numGen; g++ )
			realTimeCost += (inst.powSys->generators[g].variableCost*runParam.ED_resolution/60)*(genUsed[g][t] + overGen[g][t]);

		/* Load shedding penalty */
		for ( int d = 0; d < numBus; d++ )
			realTimeCost += loadShedPenaltyCoef*demShed[d][t];
	}
	obj = IloMinimize(env, realTimeCost);
	obj.setName(elemName);
	model.add(obj);
	realTimeCost.end();

#if defined(WRITE_PROB)
	model.exportModel("rtED.lp")
#endif

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
			for (int g = 0; g < numGen; g++) {
				for (int t = 0; t < numPeriods; t++) {
					inst.solution.g_ED[g][t0+t] = cplex.getValue(genUsed[g][t]) + cplex.getValue(overGen[g][t]);
				}
			}
		}
	}
	catch (IloException &e) {
		cout << e << endl;
	}

	return status;
}
