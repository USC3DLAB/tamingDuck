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

ofstream cplexLog;



/* Build and solve a framework for a given time series data
 *
 * Long-term UC  - deterministic
 * Short-term UC - deterministic
 * Economic dispatch - deterministic
 *
 */
int setup_DUCDED(PowSys &powSys, StocProcess &stocProc) {
	int h, t, n;

	/* logging */
	open_file(cplexLog, "cplex.log");
	
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
		
		int beginMin = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 0; mktime(timeInfo);
		
		bool status = true;

		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);

			UCmodel DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, beginMin, rep);
			
			status = DAmodel.solve();
			if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
			else		printf("Failed.\n");
			
			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				UCmodel STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, beginMin, rep);
				
				status = STmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
				else		printf("Failed.\n");


				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					
					int ED_beginPeriod = beginMin/runParam.ED_resolution;
					
					EDmodel DED(inst, ED_beginPeriod, rep);
					DED.formulate(inst, ED_beginPeriod);
					
					status = DED.solve(inst, ED_beginPeriod);
					if (status)	printf("Success (Obj= %.2f).\n", DED.getObjValue());
					else		printf("Failed.\n");
					
					/* Move to the next period */
					beginMin += runParam.baseTime;
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
					
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
	}
	cout << "------------------------------------------------------------------" << endl;

	cplexLog.close();
	
	return 0;
}

int setup_DUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath) {
	int h, t, n;

	/* logging */
	open_file(cplexLog, "cplex.log");

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
	cout << "---------- Det UC & Stoch ED Optimization / Simulation -----------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "Observation-" << rep << endl;

		/* allocate memory to hold solutions */
		Solution soln;
		soln.allocateMem(powSys.numGen, runParam.DA_horizon/runParam.baseTime);

		int beginMin = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 0; mktime(timeInfo);

		bool status;
		
		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);

			UCmodel DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, beginMin, rep);
			status = DAmodel.solve();
			if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
			else		printf("Failed.\n");

			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				UCmodel STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, beginMin, rep);
				STmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
				else		printf("Failed.\n");
				
				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);

					int ED_beginPeriod = beginMin/runParam.ED_resolution;
					
					/* Setup the ED model in cplex.concert */
					EDmodel DED(inst, ED_beginPeriod, rep);
					DED.formulate(inst, ED_beginPeriod);

					/* Translate the data structures to suit those used in 2-SD */
					status = integrateSD(inst, DED, "rted", configPath, inst.stocObserv[2], ED_beginPeriod);
					//if (status)	printf("Success (Obj= %.2f).\n", obj_val);
					if (status)	printf("Success (Obj= ??)\n");
					else		printf("Failed.\n");

					/* Move to the next period */
					beginMin += runParam.baseTime;
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
	}
	cout << "------------------------------------------------------------------" << endl;

	cplexLog.close();
	
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
		
		int beginMin = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 0; mktime(timeInfo);
		
		bool status;
		
		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
			
			SUCmaster DAmodel;
			//FIXME: Rep?
			//DAmodel.formulate(inst, DayAhead, Transmission, 0, rep);
			DAmodel.formulate(inst, DayAhead, Transmission, beginMin);
			status = DAmodel.solve();
			if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
			else		printf("Failed.\n");

			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				SUCmaster STmodel;
				//FIXME: Rep?
				//STmodel.formulate(inst, ShortTerm, Transmission, 0, rep);
				STmodel.formulate(inst, ShortTerm, Transmission, beginMin);
				status = STmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
				else		printf("Failed.\n");
				
				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					
					int ED_beginPeriod = beginMin/runParam.ED_resolution;
					
					/* Setup the ED model in cplex.concert */
					EDmodel DED(inst, ED_beginPeriod, rep);
					DED.formulate(inst, ED_beginPeriod);
					
					/* Translate the data structures to suit those used in 2-SD */
					status = integrateSD(inst, DED, "rted", configPath, inst.stocObserv[2], ED_beginPeriod);
					//if (status)	printf("Success (Obj= %.2f).\n", obj_val);
					if (status)	printf("Success (Obj= ??)\n");
					else		printf("Failed.\n");
					
					/* Move to the next period */
					beginMin += runParam.baseTime;
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
	}
	cout << "------------------------------------------------------------------" << endl;
	
	return 0;
}
