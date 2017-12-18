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

scenarios::scenarios() {}

scenarios::~scenarios() {}

scenarios::scenarios(string inputDir, string sysName) {

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
			scenType temp;
			temp = read((inputDir + sysName + "/" + rType[r] + "/" + fType[f]), ',', true, true);
			temp.name = rType[r];
			temp.type = fType[r];
			stocProc.push_back(temp);
			numStocProc++;
			cout << "Successfully read " << (inputDir + sysName + "/" + rType[r] + "/" + fType[f]).c_str() << endl;
		}
	}
}//END scenario constructor

scenType scenarios::read(string fname, char delimiter, bool readColNames, bool readRowNames) {
	ifstream fptr;
	scenType temp;
	vector<string> tokens;
	string line;
	unsigned int n;

	if ( !open_file(fptr, fname) )
		return temp;

	temp.numVars = temp.numT = 0;
	if ( readColNames ) {
		getline ( fptr, line );
		tokens = split(line, delimiter);
		for ( n = 1; n < tokens.size(); n++ )
			temp.varNames.push_back(tokens[1]);
		temp.numVars = tokens.size();
	}

	while ( getline(fptr, line) ) {
		tokens = split(line, delimiter);

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
}

/* The subroutine splits the line of type string along the delimiters into a vector of shorter strings */
vector<string> split(string &line, char delimiter) {
	stringstream ss(line);
	string item;
	vector<string> tokens;
	while (getline(ss, item, delimiter)) {
		tokens.push_back(item);
	}
	return tokens;
}
