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

#define NAMESIZE 64

class EDmodel {

public:
	EDmodel();
	~EDmodel();

	void preprocess(instance &inst, int t0);
	void formulate(instance &inst, int beginTime);
	bool solve(instance &inst, int t0);

private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	IloArray <IloNumVarArray> gen, overGen, demMet, demShed, flow, theta;

	vector<vector<double>>  busLoad;			// load at each bus and period
	vector<vector<double>>	minGenerationReq;	// minimum production requirements (obeying assumptions)
	vector<vector<double>>	maxGenerationReq;	// minimum production requirements (obeying assumptions)

	double beginTime;
};

#endif /* EDMODEL_HPP_ */
