/*
 * EDmodel.hpp
 *
 *  Created on: Dec 12, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef EDMODEL_HPP_
#define EDMODEL_HPP_

#include <ilcplex/ilocplex.h>

#include "instance.hpp"
#include "solution.hpp"

using namespace std;

//#define WRITE_PROB

#ifndef NAMESIZE
#define NAMESIZE 64
#endif

class EDmodel {

public:
	EDmodel(instance &inst, int t0, int rep);
	~EDmodel();

	void formulate(instance &inst, int t0);
	bool solve(instance &inst, int t0);
	
	double getObjValue();

	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;

	/* Rows and columns which indicate the beginning of the new stage */
	vector <string> timeCols;
	vector <string> timeRows;

	/* We assume that the randomness occurs only in the right-hand side */
	vector <string> stocRows;

private:
	IloArray <IloNumVarArray> genUsed, overGen, demMet, demShed, flow, theta, btCharge, btDischarge, btState;

	int numGen, numBus, numLine, numLoad, numPeriods, rep, numBatteries;
	vector<vector<double>> busLoad, genMin, genMax, genRampUp, genRampDown, genAvail;
	vector<double> sysLoad;
};

#endif /* EDMODEL_HPP_ */
