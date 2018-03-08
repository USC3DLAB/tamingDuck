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

#include <time.h>

#include "instance.hpp"
#include "UCmodel.hpp"
#include "EDmodel.hpp"
#include "SUC_master.hpp"
#include "solution.hpp"
#include "integrateSD.hpp"

extern runType runParam;
extern vector<Setting> settings;

ofstream timeLog;


/* Build and solve a framework for a given time series data
 *
 * Long-term UC  - deterministic
 * Short-term UC - deterministic
 * Economic dispatch - deterministic
 *
 */
int setup_DUCDED(PowSys &powSys, StocProcess &stocProc, string &RScriptsPath) {
	int h, t, n;

	double begin_t;
	
	/* logging */
	open_file(timeLog, "time.log");
	
	/* time visualization */
	time_t rawTime;
	struct tm * timeInfo;

	time ( &rawTime );
	timeInfo = localtime( &rawTime );

	string timeStamp = getCurrentDateTime();

	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, RScriptsPath);

	cout << "------------------------------------------------------------------" << endl;
	cout << "------------ Deterministic Optimization / Simulation -------------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "------------------------------------------------------------------" << endl;
		cout << "Observation-" << rep << endl;

		/* allocate memory to hold solutions */
		int beginMin = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 1; mktime(timeInfo);
		
		bool status = true;

		begin_t = get_wall_time();
		
		inst.openLogFile( ("./solver_rep" + num2str(rep) + ".log") );
		
		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
			fflush(stdout);
			
			UCmodel DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, beginMin, rep);
			
			status = DAmodel.solve();
			if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
			else		printf("Failed.\n");
			fflush(stdout);
			
			timeLog << get_wall_time() - begin_t << endl;
			
			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				fflush(stdout);
				
				UCmodel STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, beginMin, rep);
				
				status = STmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
				else		printf("Failed.\n");
				fflush(stdout);
				
				timeLog << get_wall_time() - begin_t << endl;
				
				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					fflush(stdout);
					
					int ED_beginPeriod = beginMin/runParam.ED_resolution;
					
					EDmodel DED(inst, ED_beginPeriod, rep);
					DED.formulate(inst, ED_beginPeriod);
					
					status = DED.solve(inst, ED_beginPeriod);
					if (status)	printf("Success (Obj= %.2f).\n", DED.getObjValue());
					else		printf("Failed.\n");
					fflush(stdout);
					
					timeLog << get_wall_time() - begin_t << endl;
					
					/* Move to the next period */
					beginMin += runParam.baseTime;
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
					
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
		timeLog << "------------------------------------------------------------------" << endl;
		
		inst.closeLogFile();
	}
	cout << "------------------------------------------------------------------" << endl;

	timeLog.close();
	
	return 0;
}

int setup_DUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath, string &RScriptsPath) {
	int h, t, n;

	double begin_t;
	
	/* logging */
	open_file(timeLog, "time.log");

	/* time visualization */
	time_t rawTime;
	struct tm * timeInfo;

	time ( &rawTime );
	timeInfo = localtime( &rawTime );

	string timeStamp = getCurrentDateTime();

	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, RScriptsPath);

#ifdef DEBUG_ED
	debugIntegrateSD(inst);

	/* simulate scenarios */
	inst.simulateScenarios(runParam.numTotScen, false, 0);

	printf("Running ED in debug mode.\n");
	int ED_beginPeriod = 0;
	EDmodel DED(inst, ED_beginPeriod, 0);
	DED.formulate(inst, ED_beginPeriod);

	/* translate the data structures to suit those used in 2-SD */
	double SD_objVal;
	int stat = integrateSD(inst, DED, "rted", configPath, ED_beginPeriod, SD_objVal);
	if (stat != 1)	printf("Success (Obj= %.2f).\n", SD_objVal);
	else			printf("Failed.\n");
	return 0;
#endif

	cout << "------------------------------------------------------------------" << endl;
	cout << "---------- Det UC & Stoch ED Optimization / Simulation -----------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "Observation-" << rep << endl;

		/* random sampling */
		// srand(rep);			// Commented out because: consumed only by L-shaped algorithms, which is not used in this framework.
		inst.setRSeed(rep);

		/* allocate memory to hold solutions */
		int beginMin = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 1; mktime(timeInfo);

		bool status;
		
		begin_t = get_wall_time();
		inst.openLogFile( ("./solver_rep" + num2str(rep) + ".log") );
		
		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
			fflush(stdout);
			
			/* simulate scenarios */
			inst.simulateScenarios(runParam.numTotScen, false, rep);
			
			timeLog << get_wall_time() - begin_t << endl;

			UCmodel DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, beginMin, rep);
			status = DAmodel.solve();
			if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
			else		printf("Failed.\n");
			fflush(stdout);
			
			timeLog << get_wall_time() - begin_t << endl;

			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				fflush(stdout);
				
				UCmodel STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, beginMin, rep);
				status = STmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
				else		printf("Failed.\n");
				fflush(stdout);
				
				timeLog << get_wall_time() - begin_t << endl;
				
				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					fflush(stdout);
					
					int ED_beginPeriod = beginMin/runParam.ED_resolution;
					
					/* setup the ED model in cplex.concert */
					EDmodel DED(inst, ED_beginPeriod, rep);
					DED.formulate(inst, ED_beginPeriod);

					timeLog << get_wall_time() - begin_t << endl;
					
//					/* simulate scenarios */
//					inst.simulateScenarios(runParam.numSDScen, false);
//
//					timeLog << get_wall_time() - begin_t << endl;
					
					/* translate the data structures to suit those used in 2-SD */
					double SD_objVal;
					int stat = integrateSD(inst, DED, "rted", configPath, ED_beginPeriod, SD_objVal);
					if (stat != 1)	printf("Success (Obj= %.2f).\n", SD_objVal);
					else			printf("Failed.\n");
					fflush(stdout);

					timeLog << get_wall_time() - begin_t << endl;
					
					/* move to the next period */
					beginMin += runParam.baseTime;
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
		timeLog << "------------------------------------------------------------------" << endl;
		
		inst.closeLogFile();
	}
	cout << "------------------------------------------------------------------" << endl;

	timeLog.close();
	
	return 0;
}

int setup_SUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath, string &RScriptsPath) {
	int h, t, n;
	
	double begin_t;
	
	/* logging */
	open_file(timeLog, "time.log");
	
	/* time visualization */
	time_t rawTime;
	struct tm * timeInfo;
	
	time ( &rawTime );
	timeInfo = localtime( &rawTime );
	
	string timeStamp = getCurrentDateTime();
	
	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, RScriptsPath);
	
	cout << "------------------------------------------------------------------" << endl;
	cout << "------------- Stochastic Optimization / Simulation ---------------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "Observation-" << rep << endl;
		
		/* random sampling */
		srand(rep);
		inst.setRSeed(rep);
		
		/* allocate memory to hold solutions */
		int beginMin = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 1; mktime(timeInfo);
		
		bool status;
		
		begin_t = get_wall_time();
		inst.openLogFile( ("./solver_rep" + num2str(rep) + ".log") );
		
		/* Long-term unit commitment */
		for ( h = 0; h < runParam.DA_numSolves; h++ ) {
			printf("Long-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
			fflush(stdout);
			
			/* simulate scenarios */
			inst.simulateScenarios(runParam.numTotScen, false, rep);
			
			timeLog << get_wall_time() - begin_t << endl;

			/* solve the problem */
			SUCmaster DAmodel;
			DAmodel.formulate(inst, DayAhead, Transmission, beginMin, rep);
			
			status = DAmodel.solve();
			if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
			else		printf("Failed.\n");
			fflush(stdout);

			timeLog << get_wall_time() - begin_t << endl;

			/* Short-term unit commitment */
			for ( t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				fflush(stdout);
				
//				/* simulate scenarios */
//				inst.simulateScenarios(runParam.numLSScen, false);
//
//				timeLog << get_wall_time() - begin_t << endl;
				
				/* solve the problem */
				SUCmaster STmodel;
				STmodel.formulate(inst, ShortTerm, Transmission, beginMin, rep);
				
				status = STmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
				else		printf("Failed.\n");
				fflush(stdout);
				
				timeLog << get_wall_time() - begin_t << endl;
				
				/* Economic dispatch */
				for ( n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					fflush(stdout);
					
					int ED_beginPeriod = beginMin/runParam.ED_resolution;
					
					/* setup the ED model in cplex.concert */
					EDmodel DED(inst, ED_beginPeriod, rep);
					DED.formulate(inst, ED_beginPeriod);
					
					timeLog << get_wall_time() - begin_t << endl;
					
//					/* simulate scenarios */
//					inst.simulateScenarios(runParam.numSDScen, false);
//
//					timeLog << get_wall_time() - begin_t << endl;
					
					/* translate the data structures to suit those used in 2-SD */
					double SD_objVal;
					int stat = integrateSD(inst, DED, "rted", configPath, ED_beginPeriod, SD_objVal);
					if (stat != 1)	printf("Success (Obj= %.2f).\n", SD_objVal);
					else			printf("Failed.\n");
					fflush(stdout);
					
					timeLog << get_wall_time() - begin_t << endl;
					
					/* move to the next period */
					beginMin += runParam.baseTime;
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
				}
			}
		}
		inst.printSolution( ("./" + timeStamp + "_rep" + num2str(rep)) );
		timeLog << "------------------------------------------------------------------" << endl;
		
		inst.closeLogFile();
	}
	cout << "------------------------------------------------------------------" << endl;
	
	timeLog.close();
	
	return 0;
}


int setup (PowSys &powSys, StocProcess &stocProc, string &configPath, string &RScriptsPath) {
	/* logging */
	open_file(timeLog, "time.log");
	
	/* time visualization */
	time_t rawTime;
	struct tm * timeInfo;
	time ( &rawTime );
	timeInfo = localtime( &rawTime );
	
	/* Create an instance */
	instance inst;
	inst.initialize(&powSys, &stocProc, RScriptsPath);	
	
	/* Main Body */
	switch (settings[0]) {
		case DETERMINISTIC:
			printf("%-23s%s%s\n", "Day-Ahead UC", ": ", "Deterministic");
			break;
		case STOCHASTIC:
			printf("%-23s%s%s\n", "Day-Ahead UC", ": ", "Stochastic");
			break;
		case NA:
			printf("%-23s%s%s\n", "Day-Ahead UC", ": ", "N/A");
			break;
		default:
			break;
	}
	switch (settings[1]) {
		case DETERMINISTIC:
			printf("%-23s%s%s\n", "Short-Term UC", ": ", "Deterministic");
			break;
		case STOCHASTIC:
			printf("%-23s%s%s\n", "Short-Term UC", ": ", "Stochastic");
			break;
		case NA:
			printf("%-23s%s%s\n", "Short-Term UC", ": ", "N/A");
			break;
		default:
			break;
	}
	switch (settings[2]) {
		case DETERMINISTIC:
			printf("%-23s%s%s\n", "Economic Dispatch", ": ", "Deterministic");
			break;
		case STOCHASTIC:
			printf("%-23s%s%s\n", "Economic Dispatch", ": ", "Stochastic");
			break;
		case NA:
			printf("%-23s%s%s\n", "Economic Dispatch", ": ", "N/A");
			break;
		default:
			break;
	}
	cout << "------------------------------------------------------------------" << endl;
	for (int rep = 0; rep < runParam.numRep; rep++) {
		cout << "Day: " << rep+1 << endl;
		
		/* random sampling */
		srand(rep);
		inst.setRSeed(rep);
		
		/* allocate memory to hold solutions */
		int beginMin = 0; 	timeInfo->tm_min = 0;	timeInfo->tm_hour = 1; mktime(timeInfo);
		
		bool status;
		
		double begin_t = get_wall_time();
		inst.openLogFile( ("./solver_rep" + num2str(rep+1) + ".log") );
		
		/* simulate scenarios */
		inst.simulateScenarios(runParam.numTotScen, false, rep);
		cout << endl;

		/* Long-term unit commitment */
		for (int h = 0; h < runParam.DA_numSolves; h++) {
			printf("Day-Ahead Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
			fflush(stdout);
			
			/* solve the problem */
			timeLog << get_wall_time() - begin_t << endl;
			if (settings[0] == DETERMINISTIC) {
				UCmodel DAmodel;
				DAmodel.formulate(inst, DayAhead, Transmission, beginMin, rep);
				
				status = DAmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
				else		printf("Failed.\n");
				fflush(stdout);
			}
			else if (settings[0] == STOCHASTIC) {
				SUCmaster DAmodel;
				DAmodel.formulate(inst, DayAhead, Transmission, beginMin, rep);
				
				status = DAmodel.solve();
				if (status)	printf("Success (Obj= %.2f).\n", DAmodel.getObjValue());
				else		printf("Failed.\n");
				fflush(stdout);
			}
			else {
				printf("N/A.\n");
				fflush(stdout);
			}
			timeLog << get_wall_time() - begin_t << endl;

			/* Short-term unit commitment */
			for (int t = 0; t < runParam.ST_numSolves; t++ ) {
				printf("\tShort-term Unit-Commitment (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
				fflush(stdout);
				
				/* solve the problem */
				if (settings[1] == DETERMINISTIC) {
					UCmodel STmodel;
					STmodel.formulate(inst, ShortTerm, Transmission, beginMin, rep);
					
					status = STmodel.solve();
					if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
					else		printf("Failed.\n");
					fflush(stdout);
				}
				else if (settings[1] == STOCHASTIC) {
					SUCmaster STmodel;
					STmodel.formulate(inst, ShortTerm, Transmission, beginMin, rep);
					
					status = STmodel.solve();
					if (status)	printf("Success (Obj= %.2f).\n", STmodel.getObjValue());
					else		printf("Failed.\n");
					fflush(stdout);
				}
				else {
					printf("N/A.\n");
					fflush(stdout);
				}
				timeLog << get_wall_time() - begin_t << endl;
				
				/* Economic dispatch */
				for (int n = 0; n < runParam.ED_numSolves; n++ ) {
					printf("\t\tEconomic Dispatch (%02d:%02d): ", timeInfo->tm_hour, timeInfo->tm_min);
					fflush(stdout);
					
					int ED_beginPeriod = beginMin/runParam.ED_resolution;
					
					/* setup the ED model in cplex.concert */
					EDmodel DED(inst, ED_beginPeriod, rep);
					DED.formulate(inst, ED_beginPeriod);
					
					if (settings[2] == DETERMINISTIC) {
						status = DED.solve(inst, ED_beginPeriod);
						if (status)	printf("Success (Obj= %.2f).\n", DED.getObjValue());
						else		printf("Failed.\n");
						fflush(stdout);
					}
					else if (settings[2] == STOCHASTIC) {
						/* translate the data structures to suit those used in 2-SD, execute SD */
						double SD_objVal;
						int stat = integrateSD(inst, DED, "rted", configPath, ED_beginPeriod, SD_objVal);
						if (stat != 1)	printf("Success (Obj= %.2f).\n", SD_objVal);
						else			printf("Failed.\n");
						fflush(stdout);
					}
					else {
						printf("N/A.\n");
						fflush(stdout);
					}
					timeLog << get_wall_time() - begin_t << endl;
					
					/* move to the next period */
					beginMin += runParam.baseTime;
					timeInfo->tm_min += runParam.baseTime;
					mktime(timeInfo);
				}
			}
		}
		inst.printSolution( ("./Day" + num2str(rep+1)) );
		timeLog << "------------------------------------------------------------------" << endl;
		
		inst.closeLogFile();
	}
	cout << "------------------------------------------------------------------" << endl;
	
	timeLog.close();
	
	return 0;
}
