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

#include <dirent.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
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
			temp.type = fType[f];

			sp.push_back(temp);
			numStocProc++;
			cout << "Successfully read " << (inputDir + sysName + "/" + rType[r] + "/" + fType[f]).c_str() << endl;
		}
	}

}//END scenario constructor

OneStocProc StocProcess::read(string fname, char delimiter, bool readColNames, bool readRowNames) {
	ifstream fptr;
	OneStocProc temp;
	vector<string> tokens;
	string line;
	unsigned int n;

	if ( !open_file(fptr, fname) )
		return temp;

	temp.numVars = temp.numT = 0;
	if ( readColNames ) {
		getline ( fptr, line );
		tokens = splitString(line, delimiter);
		for ( n = 1; n < tokens.size(); n++ ) {
			//cout << "First token=" << tokens[n] << " " << n << endl;
			if (tokens[n].length() > 2) tokens[n].erase(0, 1);									// remove the first \"
			//cout << "Second token=" << tokens[n] << " " << n << endl;
			if (tokens[n].length() > 2) tokens[n].erase(tokens[n].end()-1, tokens[n].end());	// remove the last \"
			//cout << "Third token=" << tokens[n] << " " << n << endl;
			//TODO: Harshaaaa
			
			temp.mapVarNamesToIndex.insert( pair<string, int> (tokens[n], n-1) );
			//cout << temp.mapVarNamesToIndex[tokens[n]] << endl;
			
			//auto it = temp.mapVarNamesToIndex.find(tokens[n]);
			//cout << it->first << " " << it->second << endl;
			
			//it = temp.mapVarNamesToIndex.find("hasan");
			//if (it == temp.mapVarNamesToIndex.end()) {
			//	cout << "Hasan is not here" << endl;
			//}
		}
		temp.numVars = (int)temp.mapVarNamesToIndex.size();
	}

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
ScenarioType createScenarioList(StocProcess *stoc, vector<int> S_indices, int T, int numVals) {
	ScenarioType observ;

	/* Compute the number of observations that can be generated */
	for ( unsigned int n = 0; n < S_indices.size(); n++ )
		if ( stoc->sp[S_indices[n]].numT/T < numVals )
			numVals = stoc->sp[S_indices[n]].numT/T;

	for ( int rep = 0; rep < numVals; rep++ ) {
		vector<vector<double>> M;
		for ( int t = 0; t < T; t++ ) {
			vector <double> vec;
			for ( unsigned int n = 0; n < S_indices.size(); n++ ) {
				for ( unsigned int m = 0; m < stoc->sp[S_indices[n]].vals[t].size(); m++)
					vec.push_back(stoc->sp[S_indices[n]].vals[t][m]);
				observ.T = T;
			}
			M.push_back(vec);
		}
		observ.vals.push_back(M);
	}

	observ.T = T;
	observ.cnt = observ.vals.size();
	observ.numOmega = observ.vals[0][0].size();
	for ( unsigned int n = 0; n < S_indices.size(); n++ ) {
		observ.name.push_back(stoc->sp[S_indices[n]].name);
		observ.mapVarNamesToIndex.insert( stoc->sp[S_indices[n]].mapVarNamesToIndex.begin(), stoc->sp[S_indices[n]].mapVarNamesToIndex.end() );
	}
	
	

	return observ;
}//END createObservList()
