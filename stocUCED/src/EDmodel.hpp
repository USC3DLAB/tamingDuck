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

	void formulate(instance &inst, int beginTime, Solution soln);
	bool solve();

private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	IloArray <IloNumVarArray> gen, dem, flow, theta;

	double beginTime;
};

#endif /* EDMODEL_HPP_ */
