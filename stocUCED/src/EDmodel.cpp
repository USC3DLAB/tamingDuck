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

void EDmodel::formulate(instance &inst, int beginTime, Solution soln) {
	int t, i, j;
	char elemName[NAMESIZE];

	/**** Decision variables *****/
	/* Generation */
	gen = IloArray< IloNumVarArray > (env, inst.powSys->numGen);
	for ( i = 0; i < inst.powSys->numGen; i++ ) {
		Generator *gptr = &(inst.powSys->generators[i]);
		sprintf(elemName, "gen[%d]", i);

		gen[i] = IloNumVarArray(env, runParam.ED_numPeriods, 0, gptr->maxCapacity, ILOFLOAT);
		gen[i].setNames(elemName); model.add(gen[i]);
	}

	/* Demand and bus angle */
	dem = IloArray< IloNumVarArray > (env, inst.powSys->numBus);
	theta = IloArray< IloNumVarArray > (env, inst.powSys->numBus);
	for ( i = 0; i < inst.powSys->numBus; i++ ) {
		Bus *bptr = &(inst.powSys->buses[i]);
		sprintf(elemName, "dem[%d]", i);
		dem[i] = IloNumVarArray(env, runParam.ED_numPeriods, 0, IloInfinity, ILOFLOAT);
		dem[i].setNames(elemName); model.add(dem[i]);

		sprintf(elemName, "theta[%d]", i);
		theta[i] = IloNumVarArray(env, runParam.ED_numPeriods, bptr->minPhaseAngle, bptr->maxPhaseAngle, ILOFLOAT);
		theta[i].setNames(elemName); model.add(theta[i]);
	}

	/* Power flows */
	flow = IloArray< IloNumVarArray > (env, inst.powSys->numLine);
	for ( i = 0; i < inst.powSys->numLine; i++ ) {
		Line *lptr = &(inst.powSys->lines[i]);
		sprintf(elemName, "flow[%d,%d]", lptr->orig->id, lptr->dest->id);

		flow[i] = IloNumVarArray(env, runParam.ED_numPeriods, lptr->minFlowLim, lptr->maxFlowLim, ILOFLOAT);
		flow[i].setNames(elemName); model.add(flow[i]);
	}

	/***** Constraints *****/
	/* Flow balance equation */
	for (t = 0; t < runParam.ED_numPeriods; t++ ) {
		for ( i = 0; i < inst.powSys->numBus; i++ ) {
			IloExpr expr (env);
			sprintf(elemName, "flowBalance[%d][%d]", t, i);

			for ( j = 0; j < inst.powSys->numGen; j++ )
				if ( inst.powSys->generators[j].connectedBus->id == i ) expr += gen[j][t];
			for ( j = 0; j < inst.powSys->numLine; j++ ) {
				if ( inst.powSys->lines[j].dest->id == i ) expr += flow[j][t];
				if ( inst.powSys->lines[j].orig->id == i ) expr -= flow[j][t];
			}
			for ( j = 0; j < inst.powSys->numBus; j++ )
				expr -= dem[j][t];

			IloConstraint c( expr == 0 ); c.setName(elemName); model.add(c);
			expr.end();
		}
	}

	/* Line power flow equation : DC approximation */
	for (t = 0; t < runParam.ED_numPeriods; t++ ) {
		for ( i = 0; i < inst.powSys->numLine; i++ ) {
			int orig = inst.powSys->lines[i].orig->id;
			int dest = inst.powSys->lines[i].dest->id;
			sprintf(elemName, "dcApprox[%d,%d][%d]", orig, dest, t);

			IloConstraint c(  flow[i][t] == inst.powSys->lines[i].susceptance*(theta[orig][t] - theta[dest][t]) ); c.setName(elemName); model.add(c);
		}
	}

	/* Generation ramping constraints */
	t = 0;
	for ( i = 0; i < inst.powSys->numGen; i++ ) {
		sprintf(elemName, "rampUp[%d][%d]", i, t);
		IloConstraint c1( gen[i][t] -  soln.g[beginTime-1][i] <= inst.powSys->generators[i].rampUpLim); c1.setName(elemName); model.add(c1);
		sprintf(elemName, "rampDown[%d][%d]", i, t);
		IloConstraint c2( inst.powSys->generators[i].rampDownLim <= gen[i][t] - soln.g[beginTime-1][i]); c2.setName(elemName); model.add(c2);
	}
	for ( t = 1; t < runParam.ED_numPeriods; t++ ) {
		for ( i = 0; i < inst.powSys->numGen; i++ ) {
			sprintf(elemName, "rampUp[%d][%d]", i, t);
			IloConstraint c1( gen[i][t] - gen[i][t-1] <= inst.powSys->generators[i].rampUpLim); c1.setName(elemName); model.add(c1);
			sprintf(elemName, "rampDown[%d][%d]", i, t);
			IloConstraint c2( inst.powSys->generators[i].rampDownLim <= gen[i][t] - gen[i][t-1]); c2.setName(elemName); model.add(c2);
		}
	}

	/***** Objective function *****/
	IloExpr dayAheadCost (env);
	IloObjective obj;
	for ( t = 0; t < runParam.ED_numPeriods; t++) {
		for ( i = 0; i < inst.powSys->numGen; i++ )
			dayAheadCost += inst.powSys->generators[i].variableCost*gen[i][t];
	}
	obj = IloMinimize(env, dayAheadCost);
	model.add(obj);
	dayAheadCost.end();

	/* Model parameters */
	cplex.extract(model);
	cplex.setParam(IloCplex::Threads, 1);

}//END formulate()
