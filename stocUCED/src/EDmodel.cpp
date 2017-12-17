///*
// * EDmodel.cpp
// *
// *  Created on: Dec 12, 2017
// *      Author: Harsha Gangammanavar
// * Institution: Southern Methodist University
// *
// * Please send your comments or bug report to harsha (at) smu (dot) edu
// *
// */
//
//#include "EDmodel.hpp"
//
//extern runType runParam;
//
//EDmodel::EDmodel() {
//	model = IloModel(env);
//	cplex = IloCplex(env);
//}
//
//EDmodel::~EDmodel() {
//	env.end();
//}
//
//void EDmodel::formulate(instance &sys1, int beginTime, solution soln) {
//	int t, i, j;
//	char elemName[NAMESIZE];
//
//	this->inst = &sys1;
//
//	/**** Decision variables *****/
//	/* Generation */
//	gen = IloArray< IloNumVarArray > (env, sys1.numGen);
//	for ( i = 0; i < sys1.numGen; i++ ) {
//		sprintf(elemName, "gen[%d]", i);
//
//		gen[i] = IloNumVarArray(env, runParam.ED_numPeriods, 0, sys1.capacity[i], ILOFLOAT);
//		gen[i].setNames(elemName); model.add(gen[i]);
//	}
//
//	/* Demand */
//	dem = IloArray< IloNumVarArray > (env, sys1.numBus);
//	for ( i = 0; i < sys1.numLoad; i++ ) {
//		sprintf(elemName, "dem[%d]", i);
//
//		//		TODO:
//		//		dem[i] = IloNumVarArray(env, runParam.ED_numPeriods, 0, demand[b], ILOFLOAT);
//		dem[i].setNames(elemName); model.add(dem[i]);
//	}
//
//	/* Power flows */
//	flow = IloArray< IloNumVarArray > (env, sys1.numLine);
//	for ( i = 0; i < sys1.numLine; i++ ) {
//		int orig = sys1.arcs[i].first;
//		int dest = sys1.arcs[i].second;
//		sprintf(elemName, "flow[%d,%d]", orig, dest);
//
//		flow[i] = IloNumVarArray(env, runParam.ED_numPeriods, sys1.min_arc_flow[i], sys1.max_arc_flow[i], ILOFLOAT);
//		flow[i].setNames(elemName); model.add(flow[i]);
//	}
//
//	/* Bus angles */
//	theta = IloArray< IloNumVarArray > (env, sys1.numBus);
//	for ( i = 0; i < sys1.numBus; i++ ) {
//		sprintf(elemName, "theta[%d]", i);
//
//		theta[i] = IloNumVarArray(env, runParam.ED_numPeriods, sys1.thetaMin[i], sys1.thetaMax[i], ILOFLOAT);
//		theta[i].setNames(elemName); model.add(theta[i]);
//	}
//
//	/***** Constraints *****/
//	/* Flow balance equation */
//	for (t = 0; t < runParam.ED_numPeriods; t++ ) {
//		for ( i = 0; i < sys1.numBus; i++ ) {
//			IloExpr expr (env);
//			sprintf(elemName, "flowBalance[%d][%d]", t, i);
//
//			for ( j = 0; j < sys1.numGen; j++ )
//				if ( sys1.generator_loc[j] == i ) expr += gen[j][t];
//			for ( j = 0; j < sys1.numLine; j++ ) {
//				if ( sys1.arcs[j].second == i ) expr += flow[j][t];
//				if ( sys1.arcs[j].first == i ) expr -= flow[j][t];
//			}
//			for ( j = 0; j < sys1.numLoad; j++ )
//				expr -= dem[j][t];
//
//			IloConstraint c( expr == 0 ); c.setName(elemName); model.add(c);
//			expr.end();
//		}
//	}
//
//	/* Line power flow equation : DC approximation */
//	for (t = 0; t < runParam.ED_numPeriods; t++ ) {
//		for ( i = 0; i < sys1.numLine; i++ ) {
//			int orig = sys1.arcs[i].first;
//			int dest = sys1.arcs[i].second;
//			sprintf(elemName, "dcApprox[%d,%d][%d]", orig, dest, t);
//
//			IloConstraint c(  flow[i][t] == sys1.susceptance[i]*(theta[orig][t] - theta[dest][t]) ); c.setName(elemName); model.add(c);
//		}
//	}
//
//	/* Generation ramping constraints */
//	t = 0;
//	for ( i = 0; i < sys1.numGen; i++ ) {
//		sprintf(elemName, "rampUp[%d][%d]", i, t);
//		IloConstraint c1( gen[i][t] -  soln.g[beginTime-1][i] <= sys1.ramp_u_lim[i]); c1.setName(elemName); model.add(c1);
//		sprintf(elemName, "rampDown[%d][%d]", i, t);
//		IloConstraint c2( sys1.ramp_d_lim[i] <= gen[i][t] - soln.g[beginTime-1][i]); c2.setName(elemName); model.add(c2);
//	}
//	for ( t = 1; t < runParam.ED_numPeriods; t++ ) {
//		for ( i = 0; i < sys1.numGen; i++ ) {
//			sprintf(elemName, "rampUp[%d][%d]", i, t);
//			IloConstraint c1( gen[i][t] - gen[i][t-1] <= sys1.ramp_u_lim[i]); c1.setName(elemName); model.add(c1);
//			sprintf(elemName, "rampDown[%d][%d]", i, t);
//			IloConstraint c2( sys1.ramp_d_lim[i] <= gen[i][t] - gen[i][t-1]); c2.setName(elemName); model.add(c2);
//		}
//	}
//
//	/***** Objective function *****/
//	IloExpr dayAheadCost (env);
//	IloObjective obj;
//	for ( t = 0; t < runParam.ED_numPeriods; t++) {
//		for ( i = 0; i < sys1.numGen; i++ )
//			dayAheadCost += sys1.var_cost[i]*gen[i][t];
//	}
//	obj = IloMinimize(env, dayAheadCost);
//	model.add(obj);
//	dayAheadCost.end();
//
//	/* Model parameters */
//	cplex.extract(model);
//	cplex.setParam(IloCplex::Threads, 1);
//
//}//END formulate()
//
//bool EDmodel::solve() {
//
//	try{
//		bool status = cplex.solve();
//		cout << "Optimization is completed with status " << cplex.getCplexStatus() << endl;
//		return status;
//	}
//	catch (IloException &e) {
//		cout << e << endl;
//		return false;
//	}
//
//
//}//END solve()
