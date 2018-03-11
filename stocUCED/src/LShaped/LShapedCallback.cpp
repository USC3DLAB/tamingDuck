//
//  LShapedCallback.cpp
//  tamingDuck
//
//  Created by Semih Atakan on 3/8/18.
//  Copyright © 2018 Semih Atakan. All rights reserved.
//

#include "LShapedCallback.hpp"
#include "SUC_master.hpp"

LShapedCallback::LShapedCallback (SUCmaster &master) {
	this->master = &master;
}

void LShapedCallback::invoke(const IloCplex::Callback::Context &context) {
	mutex.lock();	// prevent multiple threads to be below this line

//	/* Randomized Rounding Heuristic */
//	if ( context.inRelaxation() ) {
//		randomizedRounding(context);
//	}
	
	/* Benders' Algorithm */
	if ( context.inCandidate() ) {

		bool cutAdded = addLShapedCuts (context);
			
		// if no cut is added, this master-solution's recourse obj value is accurately evaluated.
		if (!cutAdded) setExpGenAmounts(context);
	}
	mutex.unlock();	// unlock, so that a new thread may process above
}

/* Something is wrong with this function
inline void LShapedCallback::randomizedRounding(const IloCplex::Callback::Context &context) {
	/* Heuristic frequency: Every 100 nodes, after node 10 *
//	if (context.getLongInfo(IloCplex::Callback::Context::Info::NodeCount) < 10
//		|| context.getLongInfo(IloCplex::Callback::Context::Info::NodeCount) % 100 != 0){
//		return;
//	}

	/* Construct a Solution *
	IloNumVarArray	vars (context.getEnv());
	IloNumArray		vals (context.getEnv());
	vector<bool> 	vals_bool;

	for (int g=0; g<master->x.getSize(); g++) {
		vars.add( master->x[g] );

		IloNumArray temp (context.getEnv());
		try {
			context.getRelaxationPoint(master->x[g], temp);
		}
		catch (IloException &e) {
			cout << e << endl;
		}
		for (int t=0; t<temp.getSize(); t++) {
			//temp[t] = round(temp[t]);
			temp[t] = (double(rand())/double(RAND_MAX) < temp[t])*1.0;
			vals_bool.push_back( round(temp[t]) );
		}
		vals.add(temp);
		temp.end();
	}

	/* Evaluate the Solution *
	try {
		// Important: Check if the produced solution is evaluated before. If it is, don't evaluate!
		// It is important that you don't visit the same solution multiple times, especially
		// if you're not evaluating them in the lazy-constraint callback
		auto it = master->evaluatedSolns.find(vals_bool);
		if (it == master->evaluatedSolns.end()) {
			master->inst->out() << "Heuristic is evaluating a new solution" << endl;
			context.postHeuristicSolution(vars, vals, 0, IloCplex::Callback::Context::SolutionStrategy::Propagate);
		} else {
			master->inst->out() << "Heuristic has re-created an observed solution" << endl;
		}
	}
	catch (IloException &e) {
		cout << e << endl;
	}
}
 */

inline void LShapedCallback::setExpGenAmounts(const IloCplex::Callback::Context &context) {
	/* record the new solution's expected initial generation amounts (to be fed to the ED problem) */
	if (master->probType == DayAhead || (master->probType == ShortTerm && master->beginMin == 0) ) {
		// initial generation amounts are only recorded in the DA-UC, and ST-UC at the planning horizon
		
		if (context.getCandidateObjective() < context.getIncumbentObjective()) {	// if this solution is improving the incumbent
			// compute the expected initial generation
			master->expInitGen = master->recourse.getExpInitGen();
		}
	}
}

inline map<vector<bool>, double>::iterator LShapedCallback::findSolution(const IloCplex::Callback::Context &context) {
	// construct the solution
	vector<bool> x (master->numGen * master->numPeriods);
	for (int g=0; g<master->numGen; g++) {
		for (int t=0; t<master->numPeriods; t++) {
			x[g*master->numPeriods + t] = context.getCandidatePoint(master->x[g][t]) > 0.5 ? 1 : 0;
		}
	}
	
	// search for the solution
	map<vector<bool>, double>::iterator it = solMap.find(x);
	
	// if not found, set its recourse obj value as infinity
	if (it == solMap.end())	{
		solMap.insert( pair<vector<bool>, double> (x, INFINITY) );
		it = solMap.find(x);
	}

	// return the iterator, so that we can chg the objective value
	return it;
}

inline bool LShapedCallback::addLShapedCuts(const IloCplex::Callback::Context &context) {
	
	/* query and record the solution */
	for (int g=0; g<master->numGen; g++) {
		context.getCandidatePoint(master->x[g], master->xvals[g]);	// retrieve generator g's solution
		for (int t=0; t<master->numPeriods; t++) {	// round the solution vector
			master->xvals[g][t] = round( master->xvals[g][t] );
		}
	}

	/* test if the solution was seen before. This could happen if:
	 * - a cut was already added, so the callback checks if there are any new cuts to add (there cannot be, since the subproblems are the same)
	 * - heuristics, warm-starts, ...
	 */
	auto solItr = findSolution(context);
	// confirm if the expectation is within relative and absolute optimality tolerances
	if (solItr->second < INFINITY) {
		double objDiff = solItr->second - context.getCandidatePoint(master->eta[0]);
		if ( objDiff/(fabs(solItr->second)+1e-14) <= RelLShapedCutOptTol || objDiff <= AbsLShapedCutOptTol ) {
			master->inst->out() << "- Solution is accepted" << endl;
			return false;
		}
	}
	
	/* set the master-problem solution in the subproblems */
	master->recourse.setMasterSoln();

	/* solve the recourse problem */
	bool isFeasible = master->recourse.solve();

	/* check if the previously-recorded objective value is improved */
	if (isFeasible && solItr->second < INFINITY) {
		/* Note: The recourse-objective may not always be accurately evaluated. The reasons are:
		 * (1) master solutions are rounded to the nearest-integer and rounded-values are fed to
		 * the subproblems,
		 * (2) numerical tolerances of the subproblems.
		 * If the objective is not improving, then there is no point of adding another cut. In
		 * fact, if you add a cut, the L-Shaped algorithm will be stuck on the same solution. */
		double objDiff = master->recourse.getObjValue() - solItr->second;
		if ( fabs(objDiff)/(fabs(solItr->second)+1e-14) <= 1e-4 || fabs(objDiff) <= 1e-6) {
			master->inst->out() << "- Cannot compute the recourse more accurately [RecObj= " << master->recourse.getObjValue() << ", PrevObj= " << solItr->second << "]. Solution is accepted" << endl;
			return false;
		}
	}
	
	/* record its objective value (for future checks) */
	if (isFeasible) {
		solItr->second = master->recourse.getObjValue();
	}

	/* check if the cut should be appended */
	double objDiff = -1;
	bool addCut;
	if (isFeasible) {
		 objDiff = master->recourse.getObjValue() - context.getCandidatePoint(master->eta[0]);
		
		// confirm if within relative and absolute optimality tolerances
		if ( objDiff/(fabs(master->recourse.getObjValue())+1e-14) > RelLShapedCutOptTol && objDiff > AbsLShapedCutOptTol ) {
			addCut = true;
		} else {
			addCut = false;
		}
	}
	else {	// infeasible
		addCut = true;
	}

	/* query and print out the recourse obj values */
	master->inst->out() << scientific << setprecision(5) << "- Obj= " << context.getCandidateObjective() << " eta= " << context.getCandidatePoint(master->eta[0]) << " RecObj= " << master->recourse.getObjValue() << " RecRelGap= " << objDiff/(fabs(master->recourse.getObjValue())+1e-14) << " RecAbsGap= " << objDiff << flush;
	
	/* prepare and add the cut */
	if (addCut) {
		// ~ Optimality Cut ~
		if (isFeasible) {
			// no multicut for now
			if (master->multicut) {
				cout << "Error: Multi-cut is not available!" << endl;
				exit(2);
			}
			
			IloNum pi_b = 0.0;
			
			// get the expected Benders' cut coefficients
			vector<vector<double>> expCutCoefs (master->numGen, vector<double> (master->numPeriods, 0.0));
			for (int s=0; s<master->rndPermutation.size(); s++) {	// for every scenario
				// get Benders' cut coefficients for scenario s
				BendersCutCoefs *cutCoefs = &(master->recourse.cutCoefs[s]);
				
				// computing the expectation, ...
				for (int g=0; g<master->numGen; g++) {
					for (int t=0; t<master->numPeriods; t++) {
						if ( fabs(cutCoefs->pi_T[g][t]) > 1e-10 ) {
							expCutCoefs[g][t] += cutCoefs->pi_T[g][t];
						}
					}
				}
				pi_b += cutCoefs->pi_b;
			}
			pi_b  *= 1.0 / (double) master->rndPermutation.size();

			// create the Benders' cut expression
			IloExpr pi_Tx (context.getEnv());
			for (int g=0; g<master->numGen; g++) {
				for (int t=0; t<master->numPeriods; t++) {
					if ( fabs(expCutCoefs[g][t]) > 1e-10 ) {
						expCutCoefs[g][t] *= 1.0 / (double) master->rndPermutation.size();
						pi_Tx += expCutCoefs[g][t] * master->x[g][t];
					}
				}
			}
			
			// add the cut
			context.rejectCandidate( master->eta[0]-pi_Tx >= pi_b );
			master->inst->out() << " (optcut)" << endl;
			
			pi_Tx.end();
		}
		else { // ~ Feasibility Cut ~
			IloExpr pi_Tx (context.getEnv());
			IloNum pi_b = 0.0;
			
			// get the infeasible scenario index
			int s = master->recourse.getInfeasScenIndex();
			
			// get Benders' cut coefficients
			BendersCutCoefs *cutCoefs = &(master->recourse.cutCoefs[s]);
			
			// create the cut
			for (int g=0; g<master->numGen; g++) {
				for (int t=0; t<master->numPeriods; t++) {
					if ( fabs(cutCoefs->pi_T[g][t] > 1e-10) ) {
						pi_Tx += cutCoefs->pi_T[g][t] * master->x[g][t];
					}
				}
			}
			pi_b += cutCoefs->pi_b;

			// add the cut
			context.rejectCandidate( -pi_Tx >= pi_b );
			master->inst->out() << " (feascut)" << endl;
			
			pi_Tx.end();
		}
	}
	else {
		master->inst->out() << " (no cut)" << endl;
	}
	
	return addCut;
}
