/*
 * setup.cpp
 *
 *  Created on: Dec 13, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "instance.hpp"
#include "UCmodel.hpp"
#include "EDmodel.hpp"
#include "solution.hpp"

extern runType runParam;

/* Build and solve a framework for a given time series data
 *
 * Long-term UC  - deterministic
 * Short-term UC - deterministic
 * Economic dispatch - deterministic
 *
 */
int setup_DUCDED(instance &inst) {
//	int h, t, n;
//
//	/* allocate memory to hold solutions */
//	Solution soln;
//	soln.allocateMem(inst.numGen, runParam.DA_horizon/runParam.ED_resolution);
//
//	/* Long-term unit commitment */
//	for ( h = 0; h < runParam.DA_numSolves; n++ ) {
//		UCmodel DAUC;
//		DAUC.formulate(inst, DayAhead, Transmission, 0);
//		DAUC.solve();
//		DAUC.printSolution();
//		DAUC.updateSoln(soln);
//
//		/* Short-term unit commitment */
//		for ( t = 0; t < runParam.ST_numSolves; t++ ) {
//			UCmodel STUC;
//			STUC.formulate(inst, ShortTerm, Transmission, n*runParam.ST_numPeriods);
//			STUC.solve();
//			STUC.printSolution();
//			STUC.updateSoln(soln);
//
//			/* Economic dispatch */
//			for ( n = 0; n < runParam.ED_numSolves; n++ ) {
//				EDmodel DED;
//				DED.formulate(inst, 0, soln);
//			}
//		}
// 	}

	return 0;
}
