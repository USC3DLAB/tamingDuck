#include "instance.hpp"
#include "./stocProcess/stoc.hpp"

extern runType runParam;

instance::instance () {}

bool instance::initialize(PowSys *powSys, StocProcess *stoc) {
	this->powSys = powSys;
	this->stoc	 = stoc;
	
	solution.allocateMem(powSys->numGen, runParam.numPeriods);

	/* Get the random realizations */
	vector<string> randName (2);
	randName[0] = "Wind";
	randName[1] = "Solar";

	vector<int> indices;
	for (auto it=stoc->mapTypeToIndex["DA"].begin(); it != stoc->mapTypeToIndex["DA"].end(); ++it) {
		if ( find(randName.begin(), randName.end(), stoc->sp[*it].name) != randName.end() ) {	// if the name is found, push it into the list of indices
			indices.push_back(*it);
		}
	}
	DA_observ = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);
	
	indices.clear();
	for (auto it=stoc->mapTypeToIndex["ST"].begin(); it != stoc->mapTypeToIndex["ST"].end(); ++it) {
		if ( find(randName.begin(), randName.end(), stoc->sp[*it].name) != randName.end() ) {	// if the name is found, push it into the list of indices
			indices.push_back(*it);
		}
	}
	ST_observ = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);
	
	indices.clear();
	for (auto it=stoc->mapTypeToIndex["RT"].begin(); it != stoc->mapTypeToIndex["RT"].end(); ++it) {
		if ( find(randName.begin(), randName.end(), stoc->sp[*it].name) != randName.end() ) {	// if the name is found, push it into the list of indices
			indices.push_back(*it);
		}
	}
	RT_observ = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);

	/* Get the deterministic realizations */
	vector<string> detName (2);
	detName[0] = "Load";

	indices.clear();
	for (auto it=stoc->mapTypeToIndex["DA"].begin(); it != stoc->mapTypeToIndex["DA"].end(); ++it) {
		if ( find(detName.begin(), detName.end(), stoc->sp[*it].name) != detName.end() ) {	// if the name is found, push it into the list of indices
			indices.push_back(*it);
		}
	}
	DA_load = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);

	indices.clear();
	for (auto it=stoc->mapTypeToIndex["ST"].begin(); it != stoc->mapTypeToIndex["ST"].end(); ++it) {
		if ( find(detName.begin(), detName.end(), stoc->sp[*it].name) != detName.end() ) {	// if the name is found, push it into the list of indices
			indices.push_back(*it);
		}
	}
	ST_load = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);

	indices.clear();
	for (auto it=stoc->mapTypeToIndex["RT"].begin(); it != stoc->mapTypeToIndex["RT"].end(); ++it) {
		if ( find(detName.begin(), detName.end(), stoc->sp[*it].name) != detName.end() ) {	// if the name is found, push it into the list of indices
			indices.push_back(*it);
		}
	}
	RT_load = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);
	
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
