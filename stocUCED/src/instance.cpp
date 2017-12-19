#include "instance.hpp"
#include "./stocProcess/stoc.hpp"

extern runType runParam;

instance::instance () {}

bool instance::initialize(PowSys *powSys, StocProcess *stoc) {
	this->powSys = powSys;
	this->stoc	 = stoc;
	
	solution.allocateMem(powSys->numGen, runParam.numPeriods);


	/* TODO: 
	 (1) Generalize the input format 
	 (2) 24 should be replaced with however time periods are needed
	 (3) Nb of scenarios
	 */
	DA_observ = createScenarioList(stoc, stoc->mapTypeToIndex["DA"], runParam.numPeriods, runParam.numRep);
	ST_observ = createScenarioList(stoc, stoc->mapTypeToIndex["ST"], runParam.numPeriods, runParam.numRep);
	RT_observ = createScenarioList(stoc, stoc->mapTypeToIndex["RT"], runParam.numPeriods, runParam.numRep);

//	bool status;
//	status = readLoadData(inputDir + sysName + "/Load/DA.csv", DA_load);
//	if (!status)	return false;
//
//	status = readLoadData(inputDir + sysName + "/Load/ST.csv", ST_load);
//	if (!status)	return false;
//
//	status = readLoadData(inputDir + sysName + "/Load/RT.csv", RT_load);
//	if (!status)	return false;
	
	return true;
}

bool instance::readLoadData(string filepath, vector<vector<double>> &load) {
	ifstream input;
	bool status = open_file(input, filepath);
	if (!status) return false;
	
	string temp_str;
	double temp_dbl;
	
	// read the headers
	safeGetline(input, temp_str);
	
	// get the # of regions
	int numRegion = (int)count(temp_str.begin(), temp_str.end(), delimiter);
	load.resize(numRegion);
	
	// read the data
	while (!input.eof()) {
		// time stamp
		getline(input, temp_str, delimiter);
		
		// regional data
		int r;
		for (r=0; r<numRegion-1; r++) {
			input >> temp_dbl;
			load[r].push_back(temp_dbl);
			move_cursor(input, delimiter);
		}
		
		// final column (separated from above to deal with eoline token)
		input >> temp_dbl;
		load[r].push_back(temp_dbl);
		safeGetline(input, temp_str);
	}
	input.close();
	
	return true;
}
