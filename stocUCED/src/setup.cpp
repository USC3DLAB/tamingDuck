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
#include "SUC_master.hpp"
#include "solution.hpp"
#include "integrateSD.hpp"

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

	/* time visualization */
	time_t rawTime;
	struct tm * timeInfo;

	time ( &rawTime );
	timeInfo = localtime( &rawTime );

	string timeStamp = getCurrentDateTime();

	/* initializations */
	// TODO: Move this initialization to a better place.
	vector<string> detElems (1);
	detElems[0] = "Load";

	vector<string> stocElems (2);
	stocElems[0] = "Solar";
	stocElems[1] = "Wind";

	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, stocElems, detElems);

	cout << "------------------------------------------------------------------" << endl;
	cout << "------------ Deterministic Optimization / Simulation -------------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "------------------------------------------------------------------" << endl;
		cout << "Observation-" << rep << endl;

		/* allocate memory to hold solutions */
		Solution soln;
		soln.allocateMem(powSys.numGen, runParam.DA_horizon/runParam.baseTime);

		int beginPeriod = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 0; mktime(timeInfo);

		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);

			UCmodel DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, 0, rep);
			DAmodel.solve();
			cout << "Success." << endl;

			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				UCmodel STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, 0, rep);
				STmodel.solve();
				cout << "Success." << endl;


				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					EDmodel DED(inst, beginPeriod, rep);
					DED.formulate(inst, beginPeriod);
					DED.solve(inst, beginPeriod);
					cout << "Success." << endl;

					/* Move to the next period */
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
					beginPeriod ++;
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
	}
	cout << "------------------------------------------------------------------" << endl;

	return 0;
}

int setup_DUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath) {
	int h, t, n;

	/* time visualization */
	time_t rawTime;
	struct tm * timeInfo;

	time ( &rawTime );
	timeInfo = localtime( &rawTime );

	string timeStamp = getCurrentDateTime();

	/* initializations */
	// TODO: Move this initialization to a better place.
	vector<string> detElems (1);
	detElems[0] = "Load";

	vector<string> stocElems (2);
	stocElems[0] = "Solar";
	stocElems[1] = "Wind";

	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, stocElems, detElems);

	cout << "------------------------------------------------------------------" << endl;
	cout << "------------ Deterministic Optimization / Simulation -------------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "Observation-" << rep << endl;

		/* allocate memory to hold solutions */
		Solution soln;
		soln.allocateMem(powSys.numGen, runParam.DA_horizon/runParam.baseTime);

		int beginPeriod = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 0; mktime(timeInfo);

		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);

			UCmodel DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, 0, rep);
			DAmodel.solve();
			cout << "Success." << endl;

			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				UCmodel STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, 0, rep);
				STmodel.solve();
				cout << "Success." << endl;
				
				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);

					/* Setup the ED model in cplex.concert */
					EDmodel DED(inst, beginPeriod, rep);
					DED.formulate(inst, beginPeriod);

					/* Translate the data structures to suit those used in 2-SD */
					integrateSD(inst, DED, "rted", configPath, inst.stocObserv[2], beginPeriod);

					cout << "Success." << endl;

					/* Move to the next period */
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
					beginPeriod ++;
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
	}
	cout << "------------------------------------------------------------------" << endl;

	return 0;
}

int setup_SUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath) {
	int h, t, n;
	
	/* time visualization */
	time_t rawTime;
	struct tm * timeInfo;
	
	time ( &rawTime );
	timeInfo = localtime( &rawTime );
	
	string timeStamp = getCurrentDateTime();
	
	/* initializations */
	// TODO: Move this initialization to a better place.
	vector<string> detElems (1);
	detElems[0] = "Load";
	
	vector<string> stocElems (2);
	stocElems[0] = "Solar";
	stocElems[1] = "Wind";
	
	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, stocElems, detElems);
	
	cout << "------------------------------------------------------------------" << endl;
	cout << "------------- Stochastic Optimization / Simulation ---------------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "Observation-" << rep << endl;
		
		/* allocate memory to hold solutions */
		Solution soln;
		soln.allocateMem(powSys.numGen, runParam.DA_horizon/runParam.baseTime);
		
		int beginPeriod = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 0; mktime(timeInfo);
		
		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
			
			SUCmaster DAmodel;
			//FIXME: Rep?
			//DAmodel.formulate(inst, DayAhead, Transmission, 0, rep);
			DAmodel.formulate(inst, DayAhead, Transmission, 0);
			DAmodel.solve();
			cout << "Success." << endl;
			
			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				SUCmaster STmodel;
				//FIXME: Rep?
				//STmodel.formulate(inst, ShortTerm, Transmission, 0, rep);
				STmodel.formulate(inst, ShortTerm, Transmission, 0);
				STmodel.solve();
				cout << "Success." << endl;
				
				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					
					/* Setup the ED model in cplex.concert */
					EDmodel DED(inst, beginPeriod, rep);
					DED.formulate(inst, beginPeriod);
					
					/* Translate the data structures to suit those used in 2-SD */
					integrateSD(inst, DED, "rted", configPath, inst.stocObserv[2], beginPeriod);
					
					cout << "Success." << endl;
					
					/* Move to the next period */
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
					beginPeriod ++;
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
	}
	cout << "------------------------------------------------------------------" << endl;
	
	return 0;
}
