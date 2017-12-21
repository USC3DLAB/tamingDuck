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
	tm nowTime;

	// TODO: Move this initialization to a better place.
	vector<string> detElems (1);
	detElems[0] = "Load";

	vector<string> stocElems (2);
	stocElems[0] = "Solar";
	stocElems[1] = "Wind";

	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, stocElems, detElems);

	cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "---------------------------------------------Deterministic optimization-simulation--------------------------------------------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "Observation-" << rep << endl;

		/* allocate memory to hold solutions */
		Solution soln;
		soln.allocateMem(powSys.numGen, runParam.DA_horizon/runParam.baseTime);

		int beginPeriod = 0; mktime(&nowTime);
		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", nowTime.tm_hour, nowTime.tm_min);
			cout << "\n\n   HARSHAAAA, what's this function above saying??" << endl;
			cout << "if you are requesting a time stamp, would this work? " << getCurrentDateTime() << endl;
			cout << "that function is in misc.cpp" << endl;
			
			UCmodel DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, 0);
			DAmodel.solve();
			cout << "Success." << endl;

			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", nowTime.tm_hour, nowTime.tm_min);
				UCmodel STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, 0);
				STmodel.solve();
				cout << "Success." << endl;

				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", nowTime.tm_hour, nowTime.tm_min);
					EDmodel DED;
					DED.formulate(inst, beginPeriod);
					DED.solve(inst, beginPeriod);
					cout << "Success." << endl;

					/* Move to the next period */
					nowTime.tm_min += runParam.baseTime;
					mktime(&nowTime);
					beginPeriod ++;
				}
			}
		}
	}
	cout << "------------------------------------------------------------------------------------------------------------------------------" << endl;

	return 0;
}
