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

#ifndef NAMESIZE
#define NAMESIZE 32
#endif

class EDmodel {

public:
	EDmodel(instance &inst, int t0, int rep);
	~EDmodel();

	void formulate(instance &inst, int t0);
	bool solve(instance &inst, int t0);

	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;

	/* Rows and columns which indicate the beginning of the new stage */
	vector <string> timeCols;
	vector <string> timeRows;

	/* We assume that the randomness occurs only in the right-hand side */
	vector <string> stocRows;

private:
	IloArray <IloNumVarArray> gen, overGen, demMet, demShed, flow, theta;

	int numGen, numBus, numLine, numLoad, numPeriods;
	vector<vector<double>> busLoad, genMin, genRampUp, genRampDown, genCap;
};

#endif /* EDMODEL_HPP_ */
