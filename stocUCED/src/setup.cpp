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
int setup_DUCDED(PowSys &powSys, StocProcess &stocProc) {
	int h, t, n;

	// TODO: Move this initialization to a better place.
	vector<string> detElems (2);
	detElems[0] = "Load";

	vector<string> stocElems (2);
	stocElems[0] = "Solar";
	stocElems[1] = "Wind";

	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, detElems, stocElems);

	//TODO: what's this cnt?
	for (int rep = 0; rep < runParam.numRep; rep++) {
		/* allocate memory to hold solutions */
		Solution soln;
		soln.allocateMem(powSys.numGen, runParam.DA_horizon/runParam.baseTime);

//		/* Long-term unit commitment */
//		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
//			UCmodel DAmodel;
//			DAmodel.formulate(inst, DayAhead, Transmission, 0);
//			DAmodel.solve();
//
//			/* Short-term unit commitment */
//			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
//				UCmodel STmodel;
//				STmodel.formulate(inst, ShortTerm, Transmission, 0);
//				STmodel.solve();
//
//				/* Economic dispatch */
//				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
//					EDmodel DED;
//					DED.formulate(inst, 0, soln);
//				}
//			}
//		}
	}

	return 0;
}
