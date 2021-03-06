/*
 * integrate.h
 *
 *  Created on: Dec 22, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef INTEGRATESD_HPP_
#define INTEGRATESD_HPP_

#include "misc.hpp"
#include "EDmodel.hpp"
#include "./sd/twoSD.h"

//#define DEBUG_ED

int integrateSD(instance &inst, EDmodel &ED, string probName, string &configPath, int t0, double &objVal);
int readConfig(string &configPath);
oneProblem *buildOneProblem(IloModel &model, IloCplex &cplex, string probName);
oneProblem *buildOneProblem_file(string probName);
timeType *buildTimeType(IloModel &model, vector<string> rowNames, vector<string> colNames);
stocType *buildStocType(IloModel &model, vector<string> stocRows, ScenarioType *stocObserv, int t0);
void debugIntegrateSD(instance &inst);

#ifdef __cplusplus
extern "C" {
void testIntegrateSD(oneProblem *orig);
}
#endif

#endif /* INTEGRATESD_HPP_ */
