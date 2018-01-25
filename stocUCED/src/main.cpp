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
#include "./powerSys/PowSys.hpp"
#include "./stocProcess/stoc.hpp"
#include <sstream>

runType runParam;

void readRunfile (string inputDir);
void parseCmdLine(int argc, const char *argv[], string &inputDir, string &configPath, string &sysName, string &setting);

int setup_DUCDED(PowSys &powSys, StocProcess &stocProc);
int setup_DUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath);
int setup_SUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath);

int main(int argc, const char * argv[]) {
	string inputDir, configPath, sysName, setting;

	/* Request for input if the default is missing */
	parseCmdLine(argc, argv, inputDir, configPath, sysName, setting);

	/* Read the configuration file */
	readRunfile (inputDir);

    /* Read the power system */
	PowSys powSys;
    powSys.readData(inputDir, sysName);

	/* checking scenario reader */
    StocProcess stocProc(inputDir, sysName);
	
	// Switch based on the chosen setting
	if ( setting == "DUC-DED" ) {
		if( setup_DUCDED(powSys, stocProc) ) {
			perror("Failed to complete the DUC-DED run.\n");
		}
	}
	else if ( setting == "DUC-SED" ) {
		if( setup_DUCSED(powSys, stocProc, configPath) ) {
			perror("Failed to complete the DUC-SED run.\n");
		}
	}
	else if ( setting == "SUC-SED" ) {
		if( setup_SUCSED(powSys, stocProc, configPath) ) {
			perror("Failed to complete the SUC-SED run.\n");
		}
	}
	else
		perror ("Unknown setting for the instance.\n");


	return 0;
}

void parseCmdLine(int argc, const char *argv[], string &inputDir, string &configPath, string &sysName, string &setting) {

	switch (argc) {
	case 2:
		inputDir = argv[1];
		cout << "Enter the path for the config.sd : ";
		cin >> configPath;
		cout << "Enter the name of the instance : ";
		cin  >> sysName;
		cout << "Enter the setting for solve (DUC-DED, DUC-SED, SUC-SED) :";
		cin  >> setting;
		break;
	case 3:
		inputDir = argv[1];
		configPath = argv[2];
		cout << "Enter the name of the instance : ";
		cin  >> sysName;
		cout << "Enter the setting for solve (DUC-DED, DUC-SED, SUC-SED) :";
		cin  >> setting;
		break;
	case 4:
		inputDir = argv[1];
		configPath = argv[2];
		sysName = argv[3];
		cout << "Enter the setting for solve (DUC-DED, DUC-SED, SUC-SED) :";
		cin  >> setting;
		break;
	case 5:
		inputDir = argv[1];
		configPath = argv[2];
		sysName = argv[3];
		setting = argv[4];
		break;
	default:
		cout << "Enter the directory : ";
		cin >> inputDir;
		cout << "Enter config path : ";
		cin >> configPath;
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
	runParam.numRep = 1;
	
	/* Read the run parameters if a run file is included in the default folder */
	if ( open_file(fptr, (inputDir + "runParameters.txt")) ) {
		while ( getline(fptr, line) ) {
			istringstream iss(line);
			if ( iss >> field1 >> field2 >> field3) {
				/* Convert all numbers into minutes */
				 if ( field3 == "hours" )
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
				 else if ( field1 == "numRep" )
					 runParam.numRep = temp;
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
	runParam.ST_numPeriods = runParam.ST_horizon/runParam.ST_resolution; runParam.ST_numSolves = runParam.DA_frequency/runParam.ST_frequency;
	runParam.ED_numPeriods = runParam.ED_horizon/runParam.ED_resolution; runParam.ED_numSolves = runParam.ST_frequency/runParam.ED_frequency;

	/* Set the base time */
	runParam.baseTime = runParam.ED_resolution;
	
	runParam.numPeriods = (int)round(runParam.horizon/runParam.baseTime);

	fptr.close();

}//END readConfig()
