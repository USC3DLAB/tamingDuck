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
string outDir;

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
    bool status = powSys.readData(inputDir, sysName);
	if (!status) {
		exit(1);
	}

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
	string tempDir;
	if (argc == 6) {
		inputDir	= argv[1];
		tempDir	    = argv[2];
		configPath	= argv[3];
		RScriptsPath= argv[4];
		sysName		= argv[5];
		setting		= argv[6];
	}
	else if (argc == 10 && strcmp(argv[6], "-setting") == 0) {
		inputDir	= argv[1];
		tempDir     = argv[2];
		configPath	= argv[3];
		RScriptsPath= argv[4];
		sysName		= argv[5];
		setting		= "custom";	// -setting
		
		settings.resize(3);
		for (int i=0; i<3; i++) {
			if ( strcmp(argv[7+i], "d") == 0 ) {
				settings[i] = DETERMINISTIC;
			} else if ( strcmp(argv[7+i], "s") == 0 ) {
				settings[i] = STOCHASTIC;
			} else if ( strcmp(argv[7+i], "na") == 0 ) {
				settings[i] = NA;
			} else {
				cout << "Wrong input" << endl;
				exit(1);
			}
		}
	}
	else {
		cout << "Missing inputs. Please provide the following in the given order:\n"
				"(1) input directory path,\n "
				"(2) output directory path,\n"
				"(3) SD config.sd path,\n"
				"(4) R scripts path,\n"
				"(5) system name,\n"
				"(6) framework setting." << endl;
		cout << "Instead of (6), you can type ""-setting"" followed by three modeling settings, i.e., ""-setting d na s (aka, deterministic N/A stochastic)""" << endl;
		exit(1);
	}
	
	/* Create output folder */
	// SA: The below commands are not properly working on Mac
	string cmdStr;
	outDir = tempDir + sysName;
	cmdStr = "mkdir " +  outDir;
	system(cmdStr.c_str());
	cout << "\nAll output files will be written to " << outDir << endl;
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
	
	runParam.resPerc_UC = 0.2;
	runParam.resPerc_ED = 0.2;
	runParam.useGenHistory = false;
	runParam.updateForecasts = false;
	
	runParam.rampingCoef = 1.0;
	runParam.renewableCoef = 1.0;
	
	runParam.storageCoef = 1.0;
	runParam.storageDev = 0.1;

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
				else if ( field1 == "resPerc_UC" )
					runParam.resPerc_UC = temp;
				else if ( field1 == "resPerc_ED" )
					runParam.resPerc_ED = temp;
				else if ( field1 == "useGenHistory" )
					runParam.useGenHistory = temp;
				else if ( field1 == "renewableCoef" )
					runParam.renewableCoef = temp;

				else if ( field1 == "storageCoef" )
					runParam.storageCoef = temp;
				else if ( field1 == "storageDev" )
					runParam.storageDev = temp;

				else if ( field1 == "rampingCoef" )
					runParam.rampingCoef = temp;
				else if ( field1 == "updateForecasts" )
					runParam.updateForecasts = temp;
				else if ( field1 == "startRep" )
					runParam.startRep = temp;
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
	if ( runParam.startRep < 1 ) {
		perror("Starting replication cannot be < 1.\n");
	}
	
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
	cout << endl;
	cout << "UC reserve percentage = " << setprecision(2) << runParam.resPerc_UC*100 << "%" << endl;
	cout << "ED reserve percentage = " << setprecision(2) << runParam.resPerc_ED*100 << "%" << endl;
	cout << "Renewable scaling coefficient = " << runParam.renewableCoef << endl;
	cout << "Ramping rate scaling coefficient = " << runParam.rampingCoef << endl;
	cout << "Storage scaling coefficient = " << runParam.storageCoef << endl;
	if (runParam.useGenHistory) cout << "Using generator histories from earlier days." << endl;
	cout << "------------------------------------------------------------------" << endl;

}//END readConfig()
