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
#include "./powerSys/PowSys.hpp"
#include "./stocProcess/stoc.hpp"
#include <sstream>

runType runParam;
vector<Setting> settings = {DETERMINISTIC, DETERMINISTIC, DETERMINISTIC};

void readRunfile (string inputDir);
void parseCmdLine(int argc, const char *argv[], string &inputDir, string &configPath, string &RScriptsPath, string &sysName, string &setting);

int setup_DUCDED(PowSys &powSys, StocProcess &stocProc, string &RScriptsPath);
int setup_DUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath, string &RScriptsPath);
int setup_SUCSED(PowSys &powSys, StocProcess &stocProc, string &configPath, string &RScriptsPath);
int setup(PowSys &powSys, StocProcess &stocProc, string &configPath, string &RScriptsPath);

int main(int argc, const char * argv[]) {
	string inputDir, configPath, RScriptsPath, sysName, setting;
	
	/* Request for input if the default is missing */
	parseCmdLine(argc, argv, inputDir, configPath, RScriptsPath, sysName, setting);
	
	/* Read the configuration file */
	readRunfile (inputDir);

    /* Read the power system */
	PowSys powSys;
    powSys.readData(inputDir, sysName);

	/* checking scenario reader */
    StocProcess stocProc(inputDir, sysName);
	
	// Switch based on the chosen setting
	if ( setting == "DUC-DED" ) {
		if( setup_DUCDED(powSys, stocProc, RScriptsPath) ) {
			perror("Failed to complete the DUC-DED run.\n");
		}
	}
	else if ( setting == "DUC-SED" ) {
		if( setup_DUCSED(powSys, stocProc, configPath, RScriptsPath) ) {
			perror("Failed to complete the DUC-SED run.\n");
		}
	}
	else if ( setting == "SUC-SED" ) {
		if( setup_SUCSED(powSys, stocProc, configPath, RScriptsPath) ) {
			perror("Failed to complete the SUC-SED run.\n");
		}
	}
	else if ( setting == "custom" ) {
		if( setup(powSys, stocProc, configPath, RScriptsPath) ) {
			perror("Failed to complete the run.\n");
		}
	}
	else {
		perror ("Unknown framework.\n");
	}

	return 0;
}

void parseCmdLine(int argc, const char *argv[], string &inputDir, string &configPath, string &RScriptsPath, string &sysName, string &setting)
{
	if (argc == 6) {
		inputDir	= argv[1];
		configPath	= argv[2];
		RScriptsPath= argv[3];
		sysName		= argv[4];
		setting		= argv[5];
	}
	else if (argc == 9 && strcmp(argv[5], "-setting") == 0) {
		inputDir	= argv[1];
		configPath	= argv[2];
		RScriptsPath= argv[3];
		sysName		= argv[4];
		setting		= "custom";	// -setting
		
		settings.resize(3);
		for (int i=0; i<3; i++) {
			if ( strcmp(argv[6+i], "deterministic") == 0 ) {
				settings[i] = DETERMINISTIC;
			} else if ( strcmp(argv[6+i], "stochastic") == 0 ) {
				settings[i] = STOCHASTIC;
			} else if ( strcmp(argv[6+i], "na") == 0 ) {
				settings[i] = NA;
			} else {
				cout << "Wrong input" << endl;
				exit(1);
			}
		}
	}
	else {
		cout << "Missing inputs. Please provide the following in the given order:\n  (1) input directory path,\n  (2) SD config.sd path,\n  (3) R scripts path,\n  (4) system name,\n  (5) framework setting." << endl;
		cout << "Instead of (5), you can type ""-setting"" followed by three modeling settings, i.e., ""-setting deterministic na stochastic""" << endl;
		exit(1);
	}
	
	/*
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
	 */

}//END parseCmdLine()

void readRunfile (string inputDir) {
	ifstream fptr;
	string	 line, field1, field2, field3;
	double	 temp;

	/* Set default values for the run parameters */
	runParam.horizon = 24*60;
	runParam.DA_horizon = 24*60;	runParam.DA_resolution = 60; runParam.DA_frequency = 24*60;
	runParam.ST_horizon = 3*60;		runParam.ST_resolution = 15; runParam.ST_frequency = 1*60;
	runParam.ED_horizon = 60;		runParam.ED_resolution = 15; runParam.ED_frequency = 15;
	runParam.numRep = 1;
	
	runParam.spinResPerc = 0.2;
	runParam.useGenHistory = false;
	
	/* Read the run parameters if a run file is included in the default folder */
	if ( open_file(fptr, (inputDir + "runParameters.txt")) ) {
		while ( getline(fptr, line) ) {
			istringstream iss(line);
			if ( iss >> field1 >> field2) {
				iss >> field3;
				/* Convert all numbers into minutes */
				if (field3 == "hours")
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
				else if ( field1 == "numTotScen" )
					runParam.numTotScen = temp;
				else if ( field1 == "numLSScen" )
					runParam.numLSScen = temp;
				else if ( field1 == "spinResPerc" )
					runParam.spinResPerc = temp;
				else if ( field1 == "useGenHistory" )
					runParam.useGenHistory = temp;
				else if ( field1 == "renewableMultiplier" )
					runParam.renewableMultiplier = temp;
				else {
					perror("Warning:: Unidentified run parameter in the file.\n");
				}
			}
		}
	}
	else
		perror("Failed to read the run parameters, using the default parameters.\n");

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
	
	/* Print configuration summary */
	cout << "------------------------------------------------------------------" << endl;
	cout << "Problem    Horizon   Resolution       Frequency" << endl;
	cout << "DA-UC     " << fixed << setprecision(0) << setw(4) << runParam.DA_horizon/60 << setw(4) << " hr" << setw(9) << runParam.DA_resolution << " min    " << "every " << setw(2) << runParam.DA_frequency/60 << setw(4) << " hr" << endl;
	cout << "ST-UC     " << fixed << setprecision(0) << setw(4) << runParam.ST_horizon/60 << setw(4) << " hr" << setw(9) << runParam.ST_resolution << " min    " << "every " << setw(2) << runParam.ST_frequency/60 << setw(4) << " hr" << endl;
	cout << "ED        " << fixed << setprecision(0) << setw(4) << runParam.ED_horizon << " min " << setw(8) << runParam.ED_resolution << " min    " << "every " << setw(2) << runParam.ED_frequency << setw(4) << " min" << endl;
	cout << endl << "Spinning reserve percentage = " << runParam.spinResPerc*100 << "%" << endl;
	cout << "Renewable-boost multiplier = " << runParam.renewableMultiplier << endl;
	if (runParam.useGenHistory) cout << "Using generator histories from earlier days." << endl;
	cout << "------------------------------------------------------------------" << endl;

}//END readConfig()
