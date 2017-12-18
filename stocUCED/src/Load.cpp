//
//  Load.cpp
//  tamingDuck
//
//  Created by Semih Atakan on 12/17/17.
//  Copyright © 2017 Semih Atakan. All rights reserved.
//

#include "Load.hpp"

Load::Load() {}

bool Load::read(string filepath) {

	ifstream input;
	bool status = open_file(input, filepath);
	if (!status) return false;
	
	string temp_str;
	double temp_dbl;
	
	// read the headers
	safeGetline(input, temp_str);
	
	// get the # of regions
	int numRegion = (int)count(temp_str.begin(), temp_str.end(), delimiter);
	regionalLoad.resize(numRegion);
	
	// read the data
	while (!input.eof()) {
		// time stamp
		getline(input, temp_str, delimiter);
		
		// regional data
		int r;
		for (r=0; r<numRegion-1; r++) {
			input >> temp_dbl;
			regionalLoad[r].push_back(temp_dbl);
			move_cursor(input, delimiter);
		}
		
		// final column (separated from above to deal with eoline token)
		input >> temp_dbl;
		regionalLoad[r].push_back(temp_dbl);
		safeGetline(input, temp_str);
	}
	input.close();

	return true;
}

