//
//  main.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "misc.hpp"
#include "config.hpp"
#include "instance.hpp"
#include "solution.hpp"
#include "UCmodel.hpp"
#include "EDmodel.hpp"
#include "PowSys.hpp"
#include "stoc.hpp"
#include <sstream>

//#include "master.hpp"

// TODO: I have changed a few *.h to *.hpp to separate the C files of SD from the cpp files here.

runType runParam;

void readRunfile (string inputDir);
void parseCmdLine(int argc, const char *argv[], string &inputDir, string &sysName, string &setting);

int setup_DUCDED(instance &inst);

int main(int argc, const char * argv[]) {
	string inputDir, sysName, setting;

	/* Request for input if the default is missing */
	parseCmdLine(argc, argv, inputDir, sysName, setting);

	/* Read the configuration file */
	readRunfile (inputDir);

    /* Read the power system */
	PowSys system;
    system.readData(inputDir, sysName);

	/* checking scenario reader */
    scenarios scen(inputDir, sysName);

    observType observ = createObservList(scen, {1,3}, 24, 10);

	/*** TEST ***/
	instance inst;
	inst.initialize(&system, &scen, inputDir, sysName);
	
	UCmodel DAmodel;
	DAmodel.formulate(inst, DayAhead, Transmission, 0);
	DAmodel.solve();
	
	UCmodel STmodel;
	STmodel.formulate(inst, ShortTerm, Transmission, 0);
	STmodel.solve();
	/*** TEST ***/
	
	// Switch based on the chosen setting
	if ( setting == "DUC-DED" ) {
		// Setup the problem instance

//		if( setup_DUCDED(inst) ) {
//			perror("Failed to complete the DUC-DED run.\n");
//		}
	}
	else if ( setting == "DUC-SED" ) {

	}
	else if ( setting == "SUC-SED" ) {

	}
	else
		perror ("Unknown setting for the instance.\n");


//	// Create a solution instance to hold solutions from deterministic
//	solution detSol;
//	detSol.allocateMem(inst.numGen, runParam.DA_horizon/runParam.ED_resolution);
//
//	/* Setting 1: Deterministic DA-UC + deterministic ST-UC + deterministic ED */
//	/* Solve the deterministic DA-UC with 24 hour horizon with a forecast. */
//	inst.read_scenarios(DayAhead);	// read day-ahead scenarios
//	UCmodel mip;
//	mip.formulate(inst, DayAhead, Transmission, 0);
//	mip.solve();
//	mip.printSolution();
//	mip.updateSoln(detSol);
//
//	/* Solve the deterministic ST-UCs with "updated" forecasts. */
//	int STUC_horizon = runParam.ST_resolution/60*runParam.ST_numPeriods;	// in hours
//	inst.read_scenarios(ShortTerm);
//	for (n = 0; n < runParam.DA_numPeriods/STUC_horizon; n++ ) {	/* N = 6 if ST-UC horizon is 4 hours */
//		UCmodel mip;
//		mip.formulate(inst, ShortTerm, Transmission, n*STUC_horizon);
//		mip.solve();
//		mip.printSolution();
//		mip.updateSoln(detSol);
//
//		for ( m = 0; m < runParam.ST_horizon/runParam.ED_horizon; m++ ) {
//			/* Option a: M = 4 if ED horizon is 1 hour with 10 minute period length and 1 hour steps. */
//			/* Option b: M = 24 if ED horizon is 1 hour with 10 minute period length and 10 minute steps. */
//			/* Solve the deterministic ED with "newest" forecast. */
//			EDmodel edMod;
//			edMod.formulate(inst, 0, detSol);
//		}
//	}

	return 0;
}

void parseCmdLine(int argc, const char *argv[], string &inputDir, string &sysName, string &setting) {

	switch (argc) {
	case 2:
		inputDir = argv[1];
		cout << "Enter the name of the instance : ";
		cin  >> sysName;
		cout << "Enter the setting for solve (DUC-DED, DUC-SED, SUC-SED) :";
		cin  >> setting;
		break;
	case 3:
		inputDir = argv[1];
		sysName = argv[2];
		cout << "Enter the setting for solve (DUC-DED, DUC-SED, SUC-SED) :";
		cin  >> setting;
		break;
	case 4:
		inputDir = argv[1];
		sysName = argv[2];
		setting = argv[3];
		break;
	default:
		cout << "Enter the directory : ";
		cin >> inputDir;
		cout << "Enter the name of the instance : ";
		cin  >> sysName;
		cout << "Enter the setting for solve (DUC-DED, DUC-SED, SUC-SED) :";
		cin  >> setting;
		break;
	}

}//END parseCmdLine()

void readRunfile (string inputDir) {
	ifstream fptr;
	string	 line, field1, field2, field3;
	double	 temp;

	/* Set default values for the run parameters */
	runParam.horizon = 24*60;
	runParam.DA_horizon = 24*60;	runParam.DA_resolution = 60; runParam.DA_frequency = 24*60;
	runParam.ST_horizon = 3*60;		runParam.ST_resolution = 15; runParam.ST_frequency = 3*60;
	runParam.ED_horizon = 60;		runParam.ED_resolution = 15; runParam.ED_frequency = 15;

	/* Read the run parameters if a run file is included in the default folder */
	if ( open_file(fptr, (inputDir + "runParameters.txt")) ) {
		while ( getline(fptr, line) ) {
			istringstream iss(line);
			if ( iss >> field1 >> field2 >> field3) {
				/* Convert all numbers into minutes */
				 if ( field3 != "minutes" )
					temp = (double) atof(field2.c_str()) * 60.0;
				else
					temp = (double) atof(field2.c_str());

				 if ( field1 == "horizon" )
					 runParam.horizon = temp;
				 else if ( field1 == "DA_time_horizon" )
					 runParam.DA_horizon = temp;
				 else if ( field1 == "DA_time_resolution" )
					 runParam.DA_resolution = temp;
				 else if ( field1 == "DA_time_frequency" )
					 runParam.DA_frequency = temp;
				 else if ( field1 == "ST_time_horizon" )
					 runParam.ST_horizon = temp;
				 else if ( field1 == "ST_time_resolution" )
					 runParam.ST_resolution = temp;
				 else if ( field1 == "ST_time_frequency" )
					 runParam.ST_frequency = temp;
				 else if ( field1 == "ED_time_horizon" )
					 runParam.ED_horizon = temp;
				 else if ( field1 == "ED_time_resolution" )
					 runParam.ED_resolution = temp;
				 else if ( field1 == "ED_time_frequency" )
					 runParam.ED_frequency = temp;
				 else {
					 perror("Warning:: Unidentified run parameter in the file.\n");
				 }
			}
		}
	}
	else
		perror("Failed to read the run parametres, using the default parameters.\n");

	/* Make sure that the time parameters make sense */
	if ( fmod(runParam.DA_horizon, runParam.DA_resolution) != 0 || fmod(runParam.DA_frequency, runParam.DA_resolution) != 0)
		perror("The long-term time parameters are not consistent.\n");
	if ( fmod(runParam.DA_horizon, runParam.ST_frequency) != 0 )
		perror("The long- and short-term time parameters are not compatible");
	if ( fmod(runParam.ST_horizon, runParam.ST_resolution) != 0 || fmod(runParam.ST_frequency, runParam.ST_resolution) != 0)
		perror("The short-term time parameters are not consistent.\n");
	if ( fmod(runParam.ST_horizon, runParam.ED_frequency) != 0 )
		perror("The short-term and economic dispatch time parameters are not compatible");
	if ( fmod(runParam.ED_horizon, runParam.ED_resolution) != 0 || fmod(runParam.ED_frequency, runParam.ED_resolution) != 0)
		perror("The economic dispatch time parameters are not consistent.\n");


	/* Compute the remaining run parameters */
	runParam.DA_numPeriods = runParam.DA_horizon/runParam.DA_resolution; runParam.DA_numSolves = runParam.horizon/runParam.DA_frequency;
	runParam.ST_numPeriods = runParam.ST_horizon/runParam.ST_resolution; runParam.ST_numSolves = runParam.DA_horizon/runParam.ST_frequency;
	runParam.ED_numPeriods = runParam.ED_horizon/runParam.ED_resolution; runParam.ED_numSolves = runParam.ST_horizon/runParam.ED_frequency;

	fptr.close();

}//END readConfig()
