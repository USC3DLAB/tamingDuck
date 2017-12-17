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

int getContents (string dir, vector<string> &contents);
int getDirs (string dir, vector<string> &subdirs);
int getFiles (string dir, vector<string> &files);

scenarios::scenarios() {}

scenarios::~scenarios() {}

scenarios::scenarios(string inputDir, string sysName) {
	vector<string> rType, fType;

	/* read all the random variable types: these are directories in the inputDir. */
	getDirs((inputDir + sysName), rType);
	if ( rType.size() == 0 )
		cout << "Warning: There are no random variables in the problem.\n";

	/* Loop through all the random variables */
	for ( unsigned int r = 0; r < rType.size(); r++ ) {

		getFiles((inputDir + sysName + "/" + rType[r]), fType);

		if (fType.size() == 0 )
			cout << "Warning: There are no forecast type data provided, setting default.\n";

		/* Loop through all the forecast types */
		for ( unsigned int f = 0; f < fType.size(); f++ ) {
			read((inputDir + sysName + "/" + rType[r] + "/" + fType[f]), ',');
		}
	}
}

bool scenarios::read(string fname, char delimiter) {
	scenType *temp;
	ifstream fptr;
	string str;

	if ( !open_file(fptr, fname) )
		return false;

	while (getline(fptr, str, delimiter))
		temp->varNames.push_back(str);

	size_t i = 0;
	while (i < temp->varNames.size() - 1) {
		if ( getline(fptr, str) ) break;

		istringstream iss( str );
		vector <string> record;

		while (iss) {
			if (!getline( iss, str, delimiter)) break;
			record.push_back(str);
		}
	}

	return true;
}//read()
