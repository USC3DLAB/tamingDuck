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

#include "EDmodel.hpp"
#include "twoSD.h"

#define TRANSLATE

int integrateSD(EDmodel ED);
oneProblem *buildOneProblem(IloModel model, IloCplex cplex, string probName);
oneProblem *buildOneProblem_file(string probName);
timeType *buildTimeType();

#ifdef __cplusplus
extern "C" {
void testIntegrateSD(oneProblem *orig);
}
#endif

#endif /* INTEGRATESD_HPP_ */
