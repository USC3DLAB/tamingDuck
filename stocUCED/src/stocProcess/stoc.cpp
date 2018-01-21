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
		getline ( fptr, line );
		tokens = splitString(line, delimiter);
		for ( n = 1; n < tokens.size(); n++ ) {
			unsigned int first = tokens[n].find_first_of('"');
			unsigned int last = tokens[n].find_last_of('"');
			if ( first != std::string::npos && last != std::string::npos ) //TODO: I'm getting an interesting warning for this function, saying it will always be true
				tokens[n] = tokens[n].substr(first+1, last-1);

			temp.mapVarNamesToIndex.insert( pair<string, int> (tokens[n], n-1) );
		}
		temp.numVars = (int) temp.mapVarNamesToIndex.size();
	}

	/* read the data */
	while ( getline(fptr, line) ) {
		tokens = splitString(line, delimiter);

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
ScenarioType createScenarioList(StocProcess *stoc, vector<int> S_indices, int lenT, int stepSize, int numVals) {
	ScenarioType observ;

	/* Compute the number of observations that can be generated */
	for ( unsigned int n = 0; n < S_indices.size(); n++ ) {
		double numValsInData = ceil(  (double) (stoc->sp[S_indices[n]].numT - (lenT-stepSize)) / (double) stepSize  );
		if (numVals > numValsInData)
			numVals = numValsInData;
	}

	for ( int rep = 0; rep < numVals; rep++ ) {
		vector<vector<double> > M;
		int offset = rep*stepSize;
		for ( int t = offset; t < offset+lenT; t++ ) {
			vector <double> vec;
			for ( unsigned int n = 0; n < S_indices.size(); n++ ) {
				for ( unsigned int m = 0; m < stoc->sp[S_indices[n]].vals[t].size(); m++)
					vec.push_back(stoc->sp[S_indices[n]].vals[t][m]);
			}
			M.push_back(vec);
		}
		observ.vals.push_back(M);
	}

	observ.T = lenT;
	observ.cnt = (int)observ.vals.size();
	observ.numVars = (int)observ.vals[0][0].size();
	
	// register the new indices (if we have solar+wind, the indices of wind generators will be numSolar+1, numSolar+2, ...)
	int offset = 0;
	for ( unsigned int n = 0; n < S_indices.size(); n++ ) {
		for (auto it = stoc->sp[S_indices[n]].mapVarNamesToIndex.begin(); it != stoc->sp[S_indices[n]].mapVarNamesToIndex.end(); ++it) {
			observ.mapVarNamesToIndex.insert( pair<string, int> (it->first, offset+it->second) );
		}
		offset += observ.mapVarNamesToIndex.size();
	}

	return observ;
}//END createObservList()
