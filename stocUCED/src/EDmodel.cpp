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
#include "./powerSys/Generator.hpp"

extern runType runParam;
extern string outDir;

EDmodel::EDmodel(instance &inst, int t0, int rep) {

	model = IloModel(env);
	cplex = IloCplex(model);

	//cplex.setOut(env.getNullStream());
	cplex.setOut( inst.out() );
	cplex.setWarning( inst.out() );


	/* Pre-process to setup some parameters: generator status (ramping and capacity), load, and stochastic generation */
	numBus = numLoad = inst.powSys->numBus;
	numGen = inst.powSys->numGen;
	numLine = inst.powSys->numLine;
	this->rep = rep;
	numPeriods = runParam.ED_numPeriods;
	numBatteries = inst.powSys->numBatteries;

	resize_matrix(busLoad, numBus, numPeriods);
	resize_matrix(genMin, numGen, numPeriods);
	resize_matrix(genMax, numGen, numPeriods);
	resize_matrix(genRampDown, numGen, numPeriods);
	resize_matrix(genRampUp, numGen, numPeriods);
	resize_matrix(genAvail, numGen, numPeriods);
	sysLoad.resize(numPeriods);

	/* Load */
	fill( sysLoad.begin(), sysLoad.end(), 0.0 );
	for (int b=0; b<numBus; b++) {
		Bus *busPtr = &(inst.powSys->buses[b]);

		// use actual values of the load (assuming that load-forecasts are perfect)
		auto it = inst.actuals.mapVarNamesToIndex.find( num2str(busPtr->regionId) );

		for (int t=0; t<numPeriods; t++) {
			busLoad[b][t] = inst.actuals.vals[rep][t0+t][it->second] * busPtr->loadPercentage;
			sysLoad[t] += busLoad[b][t];
		}
	}

	for ( int g = 0; g < numGen; g++ ) {
		Generator genPtr = inst.powSys->generators[g];

		for (int t = 0; t < numPeriods; t++ ) {
			/* If the model horizon exceeds the number of periods, then the last period commitment are repeated for the remainder of horizon */
			int idxT = min(t0+t, runParam.numPeriods-1);

			/* Make sure that the generation capacity limits are met. The minimum and maximum for generators which are turned off are set to zero during pre-processing */
			genMax[g][t] = inst.solution.x[g][idxT] * genPtr.maxCapacity;
			genMin[g][t] = inst.solution.x[g][idxT] * min(genPtr.minGenerationReq, min(genPtr.rampUpLim * runParam.ED_resolution, genPtr.rampDownLim * runParam.ED_resolution));

			// error check
			if (genMax[g][t] < 0) {
				ofstream errorlog;
				open_file(errorlog, outDir + "/error.log");
				errorlog << "Error:" << endl;
				errorlog << "genMax g: " << g << " t: " << t << " " << setprecision(10) << genMax[g][t] << endl;
				errorlog << "genMin g: " << g << " t: " << t << " " << setprecision(10) << genMin[g][t] << endl;
				errorlog << "x g: " << g << " t: " << t << " idxT " << idxT << " " << inst.solution.x[g][idxT] << endl;
				errorlog << "maxCapacity: " << genPtr.maxCapacity << endl;
				errorlog.close();
				inst.printSolution(outDir + "/infeasRTED");
			}
			if (genMin[g][t] < 0) {
				ofstream errorlog;
				open_file(errorlog, outDir + "/error.log");
				errorlog << "Error:" << endl;
				errorlog << "genMax g: " << g << " t: " << t << " " << setprecision(10) << genMax[g][t] << endl;
				errorlog << "genMin g: " << g << " t: " << t << " " << setprecision(10) << genMin[g][t] << endl;
				errorlog << "x g: " << g << " t: " << t << " idxT " << idxT << " " << inst.solution.x[g][idxT] << endl;
				errorlog << "maxCapacity: " << genPtr.maxCapacity << endl;
				errorlog.close();
				inst.printSolution(outDir + "/infeasRTED");
			}

			/* Stochastic generation set to what is available */
			auto it = inst.actuals.mapVarNamesToIndex.find(genPtr.name);
			if ( it != inst.actuals.mapVarNamesToIndex.end() ) {
				/* supply info found within the time series */

				if ( t == 0 ) {
					it = inst.actuals.mapVarNamesToIndex.find(genPtr.name);
					genAvail[g][t] = min( inst.actuals.vals[rep][t0+t][it->second], genMax[g][t] );
				}
				else {
					if (runParam.updateForecasts) {
						it = inst.meanForecast["RT"].mapVarNamesToIndex.find(genPtr.name);
						genAvail[g][t] = min(inst.meanForecast["RT"].vals[rep][t0+t][it->second], genMax[g][t]);
					}
					else {
						it = inst.meanForecast["DA"].mapVarNamesToIndex.find(genPtr.name);
						genAvail[g][t] = min(inst.meanForecast["DA"].vals[rep][t0+t][it->second], genMax[g][t]);
					}
				}
			}

			/* Ramping restrictions */
			genRampUp[g][t]		= genPtr.rampUpLim * runParam.ED_resolution;
			genRampDown[g][t]	= genPtr.rampDownLim * runParam.ED_resolution;

			// is the generator shutting down?
			if ( max(inst.solution.x[g][idxT] - inst.solution.x[g][min(idxT+1,runParam.numPeriods-1)], 0.0) > 0.5 ) {
				genRampDown[g][t] = genPtr.maxCapacity;
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

	/**** Decision variables *****
	 * IMP: Declare every variable immediately in here if the variable
	 * can appear in the first-stage problem.
	 */
	genUsed = IloArray< IloNumVarArray > (env, numGen);
	overGen = IloArray< IloNumVarArray > (env, numGen);
	demMet = IloArray< IloNumVarArray > (env, numLoad);
	demShed = IloArray <IloNumVarArray> (env, numLoad);
	theta = IloArray< IloNumVarArray > (env, numLoad);
	flow = IloArray< IloNumVarArray > (env, numLine);
	btFlow 	= IloArray<IloNumVarArray> (env, numBatteries);
	btState = IloArray<IloNumVarArray> (env, numBatteries);

	IloArray<IloNumVarArray> gamma_pos (env, numBatteries);
	IloArray<IloNumVarArray> gamma_neg (env, numBatteries);

	for ( int t = 0; t < numPeriods; t++ ) {
		/* Generation and over-generation */
		for ( int g = 0; g < numGen; g++ ) {
			if ( t == 0 ) {
				// Important: Variable declaration orders must not be altered.
				// Otherwise, SD solution-extraction needs to be updated.
				genUsed[g] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
				overGen[g] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
			}

			sprintf(elemName, "genUsed(%d)(%d)", g, t);
			genUsed[g][t].setName(elemName); model.add(genUsed[g][t]);

			if ( t == 0 && g == 0 )
				timeCols.push_back(elemName);
			else if ( t == 1 && g == 0 )
				timeCols.push_back(elemName);

			sprintf(elemName, "overGen(%d)(%d)", g, t);
			overGen[g][t].setName(elemName); model.add(overGen[g][t]);
		}

		/* Demand met, demand shed and bus angle */
		for ( int b = 0; b < numBus; b++ ) {
			Bus *bptr = &(inst.powSys->buses[b]);

			if ( t == 0 ) {
				// Important: Variable declaration orders must not be altered.
				// Otherwise, SD solution-extraction needs to be updated.
				demMet[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
				demShed[b] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
				theta[b] = IloNumVarArray(env, numPeriods, bptr->minPhaseAngle, bptr->maxPhaseAngle, ILOFLOAT);
			}

			sprintf(elemName, "demMet(%d)(%d)", b, t);
			demMet[b][t].setName(elemName); model.add(demMet[b][t]);

			sprintf(elemName, "demShed(%d)(%d)", b, t);
			demShed[b][t].setName(elemName); model.add(demShed[b][t]);

			sprintf(elemName, "theta(%d)(%d)", b, t);
			theta[b][t].setName(elemName); model.add(theta[b][t]);
		}

		/* Power flows */
		for ( int l = 0; l < numLine; l++ ) {
			Line *lptr = &(inst.powSys->lines[l]);

			if ( t == 0 ) {
				flow[l] = IloNumVarArray(env, numPeriods, lptr->minFlowLim, lptr->maxFlowLim, ILOFLOAT);
			}

			sprintf(elemName, "flow(%s)(%d)", lptr->name.c_str(), t);
			flow[l][t].setName(elemName); model.add(flow[l][t]);
		}

		/* Batteries */
		for (int bt=0; bt<numBatteries; bt++) {
			Battery *btPtr = &(inst.powSys->batteries[bt]);
			if ( t == 0 ) {
				btFlow[bt] 	= IloNumVarArray(env, numPeriods, -btPtr->maxCapacity, btPtr->maxCapacity, ILOFLOAT);
				btState[bt] = IloNumVarArray(env, numPeriods, 0, btPtr->maxCapacity, ILOFLOAT);
				gamma_pos[bt] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
				gamma_neg[bt] = IloNumVarArray(env, numPeriods, 0, IloInfinity, ILOFLOAT);
			}

			sprintf(elemName, "v(%d)(%d)", bt, t);
			btFlow[bt][t].setName(elemName);
			model.add(btFlow[bt][t]);

			sprintf(elemName, "I(%d)(%d)", bt, t);
			btState[bt][t].setName(elemName);
			model.add(btState[bt][t]);

			sprintf(elemName, "gamma_pos(%d)(%d)", bt, t);
			gamma_pos[bt][t].setName(elemName);
			model.add(gamma_pos[bt][t]);

			sprintf(elemName, "gamma_neg(%d)(%d)", bt, t);
			gamma_neg[bt][t].setName(elemName);
			model.add(gamma_neg[bt][t]);
		}
	}

	/***** Constraints *****
	 * IMP: Declare every constraint immediately in here, if the
	 * constraint may appear in the first-stage.
	 */
	for (int t = 0; t < numPeriods; t++ ) {
		/* Flow balance equation */
		for ( int b = 0; b < numBus; b++ ) {
			IloExpr expr (env);
			sprintf(elemName, "flowBalance(%d)(%d)", b, t);

			// production
			for ( int g = 0; g < numGen; g++ )
				if ( inst.powSys->generators[g].connectedBus->id == b ) expr += genUsed[g][t];

			// storage
			Bus* busPtr = &(inst.powSys->buses[b]);
			for (int bt=0; bt < (int) busPtr->connectedBatteries.size(); bt++) {
				expr -= btFlow[ busPtr->connectedBatteries[bt]->id ][t];
			}

			// in/out flow
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
			sprintf(elemName, "dcApprox(%s)(%d)", inst.powSys->lines[l].name.c_str(), t);

			IloExpr expr (env);
			expr = flow[l][t] - inst.powSys->lines[l].susceptance*(theta[orig][t] - theta[dest][t]);

			IloConstraint c(expr == 0); c.setName(elemName); model.add(c);
		}

		/* Generation ramping constraints: applicable only if the generator continues to be ON, i.e., x[g][t] = 1. */
		int idx = min(t0+t, runParam.numPeriods-1);

		if ( t == 0 ) {
			// Note: Generation amounts at time t-1 comes from:
			// a) If t-1 is earlier than the planning horizon, 1st-period generation levels of the DA-UC problem,
			// b) Otherwise, the ED solutions obtained earlier.
			double prevGen;

			for (int g=0; g<numGen; g++) {
				if (inst.solution.x[g][idx] < 0.5)	continue;		// generator is not operational

				// determine previous generation level
				if ( t0 == 0 ) {
					if (runParam.useGenHistory && inst.solList.size() > 0) {
						prevGen = inst.solList.back().g_ED[g][runParam.numPeriods-1];
					} else {
						prevGen = inst.solution.g_UC[g][t0] * 1.0;	// all generators are assumed to be operational
					}
				} else {
					prevGen = inst.solution.g_ED[g][t0-1] * round(inst.solution.x[g][t0-1]);	// the latter is to prevent numerical errors
				}

				// add the constraints
				/* ramp-up */
				sprintf(elemName, "rampUp(%d)(%d)", g, t);
				IloConstraint c1( genUsed[g][t] + overGen[g][t] - prevGen <= genRampUp[g][t]);
				c1.setName(elemName); model.add(c1);

				/* ramp-down */
				sprintf(elemName, "rampDown(%d)(%d)", g, t);
				IloConstraint c2( prevGen - genUsed[g][t] - overGen[g][t] <= genRampDown[g][t] );	// the latter
				c2.setName(elemName); model.add(c2);
			}
		}
		else	/* The remaining periods of the ED horizon */
		{
			for ( int g = 0; g < numGen; g++ ) {
				if (inst.solution.x[g][idx] < 0.5)	continue;	// generator is not operational

				/* ramp-up */
				sprintf(elemName, "rampUp(%d)(%d)", g, t);
				IloConstraint c1( genUsed[g][t] + overGen[g][t] - genUsed[g][t-1] - overGen[g][t-1] <= genRampUp[g][t]);
				c1.setName(elemName); model.add(c1);

				/* ramp-down */
				sprintf(elemName, "rampDown(%d)(%d)", g, t);
				IloConstraint c2( genUsed[g][t-1] + overGen[g][t-1] - genUsed[g][t] - overGen[g][t] <= genRampDown[g][t]);
				c2.setName(elemName); model.add(c2);
			}
		}

		for ( int g = 0; g < numGen; g++ ) {
			Generator genPtr = inst.powSys->generators[g];

			/* Stochastic generation set to what is available */
			auto it = inst.actuals.mapVarNamesToIndex.find(genPtr.name);
			bool stocSupply = (it != inst.actuals.mapVarNamesToIndex.end()) ? true : false;

			if ( !stocSupply ) {
				/* Deterministic-supply generator */
				if (genPtr.isMustUse) {
					sprintf(elemName, "detAvail(%d)(%d)", g, t);
					IloConstraint c (genUsed[g][t] == genMax[g][t]); c.setName(elemName); model.add(c);
					overGen[g][t].setUB(0);
				}
				else {
					sprintf(elemName, "maxGen(%d)(%d)", g, t);
					IloConstraint c1( genUsed[g][t] + overGen[g][t] <= genMax[g][t]);
					c1.setName(elemName); model.add(c1);

					if (genMin[g][t] >= EPSzero) {	// if no min-generation requirement, skip the constraint
						sprintf(elemName, "minGen(%d)(%d)", g, t);
						IloConstraint c2( genUsed[g][t] + overGen[g][t] >= genMin[g][t]);
						c2.setName(elemName); model.add(c2);
					}
				}
			}
			else {
				/* Stochastic-supply generator */
				sprintf(elemName, "stocAvail(%d)(%d)", g, t);
				if (genPtr.isMustUse) {
					IloConstraint c (genUsed[g][t] == genAvail[g][t]); c.setName(elemName); model.add(c);
					overGen[g][t].setUB(0);
				} else {
					IloConstraint c (genUsed[g][t] + overGen[g][t] == genAvail[g][t]); c.setName(elemName); model.add(c);
				}

				if ( t != 0 )
					stocRows.push_back(elemName);
			}
		}

		/* Demand consistency */
		for ( int d = 0; d < numBus; d++ ) {
			sprintf(elemName, "demConsist(%d)(%d)", d, t);
			if (t==0) {
				IloConstraint c(demMet[d][t] + demShed[d][t] == busLoad[d][t]); c.setName(elemName); model.add(c);
			} else {
				/* Node-wise spinning reserves */
				IloConstraint c(demMet[d][t] + demShed[d][t] == busLoad[d][t]*(1+runParam.resPerc_ED)); c.setName(elemName); model.add(c);
			}
		}

		// battery state
		for (int bt = 0; bt < numBatteries; bt++) {
			Battery* batPtr = &inst.powSys->batteries[bt];

			double dissipationCoef = pow(batPtr->dissipationCoef, runParam.ED_resolution/60.0);
			double conversionLossCoef = pow(batPtr->conversionLossCoef, runParam.ED_resolution/60.0);

			if (t == 0) {
				// determine previous battery state
				double initBtState = 0;
				if (t0 == 0) {
					// TODO: What are the two cases here?
					if (runParam.useGenHistory && inst.solList.size() > 0) {
						initBtState = inst.solList.back().btState_ED[bt][runParam.numPeriods-1];
					} else {
						initBtState = inst.powSys->batteries[bt].maxCapacity / 2.0;
					}
				} else {
					initBtState = inst.solution.btState_ED[bt][t0-1];
				}
				IloConstraint c (0 == -btState[bt][t] + initBtState * dissipationCoef + btFlow[bt][t] * conversionLossCoef);
				sprintf(elemName, "Bt_%d_%d", bt, t); c.setName(elemName); model.add(c);
			}
			else {
				IloConstraint c (0 == -btState[bt][t] + btState[bt][t-1] * dissipationCoef + btFlow[bt][t] * conversionLossCoef);
				sprintf(elemName, "Bt_%d_%d", bt, t); c.setName(elemName); model.add(c);
			}

			/* TODO: battery state deviation from UC targets */
			double target = 0;
			int index;
			if ( t0+t > (1+t0/runParam.ST_numPeriods) * runParam.ST_numPeriods - 1) {
				index = (1+t0/runParam.ST_numPeriods) * runParam.ST_numPeriods - 1;
			} else {
				index = t0+t;
			}

			index = min(index, runParam.numPeriods-1);
			target = inst.solution.btState_UC[bt][index];

			IloConstraint c( btState[bt][t] + gamma_neg[bt][t] - gamma_pos[bt][t] == target );
			sprintf(elemName, "BtDev(%d)(%d)", bt, t);
			c.setName(elemName);
			model.add(c);
		}
	}

	/** Rampability to UC generation levels **/
	/* IMP: These variables and constraints only appear in the second-stage
	 * problem, therefore they can be declared after the above loop.
	 */
	IloNumVarArray delta_pos (env, numGen, 0, IloInfinity, ILOFLOAT);	// positive deviations from settled DA-UC generation amounts
	IloNumVarArray delta_neg (env, numGen, 0, IloInfinity, ILOFLOAT);	// negative deviations from settled DA-UC generation amounts	
	for (int g=0; g<numGen; g++) {
		sprintf(elemName, "delta_pos(%d)", g);
		delta_pos[g].setName(elemName);

		sprintf(elemName, "delta_neg(%d)", g);
		delta_neg[g].setName(elemName);
	}
	model.add(delta_pos); model.add(delta_neg);

	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);

		if (genPtr->type != Generator::SOLAR && genPtr->type != Generator::WIND) {
			int t = numPeriods-1;
			int tprime = numPeriods;

			double target = 0;
			int index;
			if (t0+tprime > (1+t0/runParam.ST_numPeriods) * runParam.ST_numPeriods - 1) {
				index = (1+t0/runParam.ST_numPeriods) * runParam.ST_numPeriods - 1;
			} else {
				index = t0+tprime;
			}

			index = min(index, runParam.numPeriods-1);
			target = inst.solution.g_UC[g][index];

			/* ramp-up */
			IloConstraint c1( target - genUsed[g][t] - overGen[g][t] - delta_pos[g] <= genPtr->rampUpLim * runParam.ED_resolution);
			sprintf(elemName, "UPrampability(%d)(%d)", g, t);
			c1.setName(elemName);
			model.add(c1);

			/* ramp-down */
			IloConstraint c2( genUsed[g][t] + overGen[g][t] - target - delta_neg[g] <= genPtr->rampDownLim * runParam.ED_resolution);
			sprintf(elemName, "DNrampability(%d)(%d)", g, t);
			c2.setName(elemName);
			model.add(c2);
		}
	}

	/***** Objective function *****/
	IloExpr realTimeCost (env);
	IloObjective obj;
	sprintf(elemName, "realTimeCost");
	for ( int t = 0; t < numPeriods; t++) {
		/* Generation cost */
		for ( int g = 0; g < numGen; g++ ) {
			realTimeCost += (inst.powSys->generators[g].variableCost*runParam.ED_resolution/60.0)*(genUsed[g][t] + overGen[g][t]);
			if ( inst.powSys->generators[g].type == Generator::SOLAR || inst.powSys->generators[g].type == Generator::WIND ) {
				realTimeCost += renCurtailPenaltyCoef * overGen[g][t];
			}
			else {
				realTimeCost += overGenPenaltyCoef*overGen[g][t];
			}
		}

		/* Load shedding penalty */
		for ( int d = 0; d < numBus; d++ )
			realTimeCost += loadShedPenaltyCoef*demShed[d][t];
	}

	/* Deviation penalties */
	for (int g=0; g<numGen; g++) {
		if ( inst.powSys->generators[g].type != Generator::SOLAR && inst.powSys->generators[g].type != Generator::WIND ) {
			realTimeCost += 1000*(delta_pos[g] + delta_neg[g]);
		}
	}
	for (int bt = 0; bt < numBatteries; bt++) {
		for (int t = 0; t < numPeriods; t++) {
			realTimeCost += (overGenPenaltyCoef+renCurtailPenaltyCoef)/2.0 * runParam.storageCoef * (gamma_pos[bt][t] + gamma_neg[bt][t]);
		}
	}

	obj = IloMinimize(env, realTimeCost);
	obj.setName(elemName);
	model.add(obj);
	realTimeCost.end();

#if defined(WRITE_PROB)
	cplex.exportModel("rtED.lp");
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

		if (status) {
			double totLoadShed = 0;
			for (int b=0; b<numBus; b++) {
				for (int t=0; t<1; t++) {
					try {
						if ( cplex.getValue(demShed[b][t]) > 1e-6 ) {
							totLoadShed += cplex.getValue(demShed[b][t]);
						}
					}
					catch (IloException &e) { }
				}
			}
			if (totLoadShed > 0) {
				cout << "[LS! " << totLoadShed << " MWs] ";
			}
		}

		if (!status) {
			string fname = outDir + "infeasible_RTED.lp";
			cplex.exportModel(fname.c_str());
			exit(5);
		}

		// record the solution
		if (status) {
			for (int t = 0; t < numPeriods; t++) {
				if ( t0 + t < runParam.numPeriods ) {	// do not record for periods that exceed the planning horizon
					// used and over-generation amounts
					for (int g = 0; g < numGen; g++) {
						inst.solution.usedGen_ED[g][t0+t] = cplex.getValue(genUsed[g][t]);
						inst.solution.overGen_ED[g][t0+t] = cplex.getValue(overGen[g][t]);
						inst.solution.g_ED[g][t0+t] = inst.solution.usedGen_ED[g][t0+t] + inst.solution.overGen_ED[g][t0+t];
						inst.solution.g_ED[g][t0+t] = max(0.0, inst.solution.g_ED[g][t0+t]);	// numerical corrections
					}

					// load-shedding amounts
					for (int b = 0; b < numBus; b++) {
						inst.solution.loadShed_ED[b][t0+t] = cplex.getValue(demShed[b][t]);
					}

					// storage
					for (int bt=0; bt<numBatteries; bt++) {
						inst.solution.btState_ED[bt][t0+t] = cplex.getValue(btState[bt][t]);
						inst.solution.btFlow_ED[bt][t0+t] = cplex.getValue(btFlow[bt][t]);
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

/****************************************************************************
 * getObjValue
 * - Returns the objective value of the last solve.
 ****************************************************************************/
double EDmodel::getObjValue() {
	return cplex.getObjValue();
}
