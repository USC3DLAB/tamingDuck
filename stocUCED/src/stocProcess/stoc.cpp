/*
 * stoc.cpp
 *
 *  Created on: Dec 14, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "../misc.hpp"
#include "stoc.hpp"

extern runType runParam;

StocProcess::StocProcess() {}

StocProcess::~StocProcess() {}

StocProcess::StocProcess(string inputDir, string sysName) {

	numStocProc = 0;

	/* read all the random variable types: these are directories in the inputDir. */
	vector<string> rType;
	getDirs((inputDir + sysName), rType);
	if ( rType.size() == 0 )
		cout << "Warning: There are no random variables in the problem.\n";

	/* Loop through all the random variables */
	for ( unsigned int r = 0; r < rType.size(); r++ ) {
		if (rType[r] == "stocProcesses") {
			continue;
		}
		
		vector<string> fType;

		getFiles((inputDir + sysName + "/" + rType[r]), fType);

		if (fType.size() == 0 )
			cout << "Warning: There is no type data provided, setting default.\n";

		/* Loop through all the forecast types */
		for ( unsigned int f = 0; f < fType.size(); f++ ) {
			OneStocProc temp;
			temp = read((inputDir + sysName + "/" + rType[r] + "/" + fType[f]), ',', true, true);
			temp.name = rType[r];
			temp.type = fType[f].substr(0, fType[f].find_last_of("."));
			
			sp.push_back(temp);
			
			// check if the key exists, create a new key, if it doesn't exist
			auto it = mapTypeToIndex.find(temp.type);
			if (it == mapTypeToIndex.end()) { mapTypeToIndex.insert( pair<string, vector<int> > (temp.type, vector<int> ())); }

			vector<int> *vecPtr = &(mapTypeToIndex[temp.type]);
			vecPtr->push_back(numStocProc);
			
			numStocProc++;

			cout << "Successfully read " << (inputDir + sysName + "/" + rType[r] + "/" + fType[f]).c_str() << endl;
		}
	}

}//END scenario constructor

OneStocProc StocProcess::read(string fname, char delimiter, bool readColNames, bool readRowNames) {
	ifstream fptr;
	OneStocProc temp;
	vector<string> tokens;
	string line, str;
	unsigned int n;

	/* open the file */
	if ( !open_file(fptr, fname) )
		return temp;

	/* read the column names if they exist */
	temp.numVars = temp.numT = 0;
	if ( readColNames ) {
		safeGetline(fptr, line);
		tokens = splitString(line, delimiter);
		for ( n = 1; n < tokens.size(); n++ ) {
			size_t first = tokens[n].find_first_of('"');
			size_t last = tokens[n].find_last_of('"');
			if ( first != std::string::npos && last != std::string::npos )
				tokens[n] = tokens[n].substr(first+1, last-1);

			temp.mapVarNamesToIndex.insert( pair<string, int> (tokens[n], n-1) );
		}
		temp.numVars = (int) temp.mapVarNamesToIndex.size();
	}

	/* read the data */
	while ( safeGetline(fptr, line) ) {
		tokens = splitString(line, delimiter);
		
		if (tokens.size() == 0)	break;	// eof

		n = 0;
		if ( readRowNames )
			temp.rowNames.push_back(tokens[n++]);


		{
			vector<double> vec;
			for ( n = 1; n < tokens.size(); n++ )
				vec.push_back(atof(tokens[n].c_str()));
			temp.vals.push_back(vec);
		}

		temp.numT++;
	}

	return temp;
}//END scenarios::read

/* This subroutine creates a list of _numVals_ observations for stochastic processes indexed by S_indices of time duration _T_ and returns a observType structure output. */
ScenarioType createScenarioList(StocProcess *stoc, vector<int> S_indices, int lenT, int stepSize, int &numVals, int dataPeriodLengthInMins) {
	ScenarioType observ;

	int numSubHours = dataPeriodLengthInMins / runParam.baseTime;
	lenT		*= 60.0/dataPeriodLengthInMins;
	stepSize	*= 60.0/dataPeriodLengthInMins;
	
	/* Compute the number of observations that can be generated */
	for ( unsigned int n = 0; n < S_indices.size(); n++ ) {
		double numValsInData = floor(  (double) (stoc->sp[S_indices[n]].numT - (lenT-stepSize)) / (double) stepSize  );
		if (numVals > numValsInData)
			numVals = numValsInData;
	}
	
	for (int rep=0; rep<numVals; rep++) {
		vector<vector<double>> M;
		int offset = rep*stepSize;
		
		for (int h=offset; h<offset+lenT; h++) {
			// Important: I'm assuming hourly data in here
			for (int subhour = 0; subhour<numSubHours; subhour++) {
				vector<double> vec;
				for (unsigned int n=0; n<S_indices.size(); n++) {
					for (unsigned int m=0; m<stoc->sp[S_indices[n]].vals[h].size(); m++) {
						vec.push_back(stoc->sp[S_indices[n]].vals[h][m]);
					}
				}
				M.push_back(vec);
			}
		}
		observ.vals.push_back(M);
	}

	observ.T = (int)observ.vals[0].size();
	observ.cnt = (int)observ.vals.size();
	observ.numVars = (int)observ.vals[0][0].size();
	
	// register the new indices (if we have solar+wind, the indices of wind generators will be numSolar+1, numSolar+2, ...)
	int offset = 0;
	for ( unsigned int n = 0; n < S_indices.size(); n++ ) {
		for (auto it = stoc->sp[S_indices[n]].mapVarNamesToIndex.begin(); it != stoc->sp[S_indices[n]].mapVarNamesToIndex.end(); ++it) {
			observ.mapVarNamesToIndex.insert( pair<string, int> (it->first, offset+it->second) );
		}
		offset += stoc->sp[S_indices[n]].numVars;
	}

	return observ;
}//END createObservList()

#ifdef _SAMPLE_USING_R
/****************************************************************************
 * createScenarioList
 * - Simulates stochastic data-elements of the formulations using the 
 provided R scripts, records them into a ScenarioType
 ****************************************************************************/
ScenarioType createScenarioList(RInside &R, bool fitModel, string &dataFolder, vector<string> &stocElems, int simLengthDays, int &numScen, int rep) {
	ScenarioType observ;
	
	R.parseEval("dataTypes <- NULL");				// set the data types
	for (int i=0; i<stocElems.size(); i++) {
		char buffer[100];
		sprintf(buffer, "dataTypes[%d] <- \"%s\"", i+1, stocElems[i].c_str());
		R.parseEval(buffer);
	}
	
	R["fileName"]		= "DA.csv";						// set the file name
	R["fitModel"]		= fitModel ? "TRUE" : "FALSE";	// fit the model?
	R["numScenarios"]	= numScen;
	R["simLength"]		= simLengthDays;
	
	if ( numScen <= 5000 ) {
		
		cout << endl << "**** Reading " << numScen << " scenarios from .csv files ****" << endl;
		
		// read from memory
		observ.vals.resize( numScen );
		for (int s=0; s<numScen; s++) {
			resize_matrix(observ.vals[s], 108, 92);
		}
			
		char fname[256];
		ifstream input;
		for (int s=0; s<numScen; s++) {
			sprintf(fname, "./Simulations/sim_rep%d_scen%d.csv", rep+1, s+1);
			open_file(input, fname);
			
			// skip the header
			string temp_str;
			safeGetline(input, temp_str);
			
			// read the data
			for (int t=0; t<108; t++) {
				
				// skip the row header
				move_cursor(input, delimiter);
				
				for (int g=0; g<92; g++) {
					input >> observ.vals[s][t][g];
					if (g != 91)	move_cursor(input, delimiter);
					else			move_cursor(input, '\n');
				}
			}
			input.close();
		}
	}
	else {

		cout << endl << "**** ERROR: NOT READING ANY SCENARIOS ****" << endl;
		
	R.parseEval("source(paste(RScriptsPath, \"runScript.R\", sep=\"\"))");				// execute the script
	
	Rcpp::NumericVector scenList = R.parseEval("scenarios");
	int numComponents	= R.parseEval("numComponents");
	int nbPeriods		= R.parseEval("simLength * simFrequency");
	
	// resize 3D-array
	observ.vals.resize( numScen );
	for (int s=0; s<numScen; s++) {
		resize_matrix(observ.vals[s], nbPeriods, numComponents);
	}
	
	// read into observ.vals, make sure the values do not exceed capacity
	int i=0;
	for (int s=0; s<numScen; s++) {
		for (int g=0; g<numComponents; g++) {
			for (int t=0; t<nbPeriods; t++) {
				observ.vals[s][t][g] = scenList[i++];
			}
		}
	}
	
	}
	
	observ.T = (int)observ.vals[0].size();
	observ.cnt = (int)observ.vals.size();
	observ.numVars = (int)observ.vals[0][0].size();
	
	// register the indices
	Rcpp::StringVector colNames = R.parseEval("colNames");
	for (unsigned int n=0; n<colNames.size(); n++) {
		observ.mapVarNamesToIndex.insert( pair<string, int> ((string)(colNames[n]), n) );
	}
	
	return observ;
}
#endif

/****************************************************************************
 * createScenarioList
 * - Reads simulated stochastic data-elements of the formulations using the
 provided R scripts, records them into a ScenarioType
 * - Be careful of the order of stocElems vector!
 ****************************************************************************/
ScenarioType createScenarioList(string &simulationsFolder, vector<string> &stocElems, int lenT, int &numScen, int rep) {
	ScenarioType observ;
	
	cout << endl << "**** Reading " << numScen << " scenarios from .csv files ****" << endl;
	cout << "**** Warning: stocElems must be provided in correct order. Num of stocElems are hardcoded. stocElems' names are hardcoded." << endl;
	
	int nGens = 31;
	int nPeriods = 108;
	
	// read from memory
	observ.vals.resize( numScen );
	for (int s=0; s<numScen; s++) {
		resize_matrix(observ.vals[s], nPeriods, nGens);
	}
	
	char fname[256];
	ifstream input;
	for (int s=0; s<numScen; s++) {
		sprintf(fname, "%ssim_rep%d_scen%d.csv", simulationsFolder.c_str(), rep+1, s+1);
		open_file(input, fname);
		
		// skip the header
		string temp_str;
		safeGetline(input, temp_str);
		
		// read the data
		for (int t=0; t<nPeriods; t++) {
			
			// skip the row header
			move_cursor(input, delimiter);
			
			for (int g=0; g<nGens; g++) {
				input >> observ.vals[s][t][g];
				if (g != nGens-1)	move_cursor(input, delimiter);
				else				move_cursor(input, '\n');
			}
		}
		input.close();
	}

	observ.T = (int)observ.vals[0].size();
	observ.cnt = (int)observ.vals.size();
	observ.numVars = (int)observ.vals[0][0].size();
	
	
	// register the indices
	vector<string> colNames;
	colNames.push_back("Solar-bus015");
	colNames.push_back("Solar-bus018");
	colNames.push_back("Solar-bus019");
	colNames.push_back("Solar-bus032");
	colNames.push_back("Solar-bus054");
	colNames.push_back("Solar-bus055");
	colNames.push_back("Solar-bus069");
	colNames.push_back("Solar-bus076");
	colNames.push_back("Solar-bus092");
	colNames.push_back("Solar-bus100");
	colNames.push_back("Solar-bus104");
	colNames.push_back("Solar-bus105");
	colNames.push_back("Solar-bus110");
	colNames.push_back("Solar-bus112");
	for (int i=0; i<17; i++) {
		char genName[256];
		sprintf(genName, "Wind%02d", i+1);
		colNames.push_back(genName);
	}

	for (unsigned int n=0; n<colNames.size(); n++) {
		observ.mapVarNamesToIndex.insert( pair<string, int> ((string)(colNames[n]), n) );
	}
	cout << "**** Reading is completed ****" << endl;
	
	return observ;
}


