#include "instance.hpp"
#include "misc.hpp"

// TODO: SMH: does this need to be in here?
extern runType runParam;

instance::instance () {}

bool instance::initialize(PowSys *powSys, scenarios *stoc, string inputDir, string sysName) {
	this->powSys = powSys;
	this->stoc	 = stoc;
	
	solution.allocateMem(powSys->numGen, (int)round(runParam.horizon/runParam.ED_resolution));

	bool status;
	status = readLoadData(inputDir + sysName + "/Load/DA.csv", DA_load);
	if (!status)	return false;
	
	status = readLoadData(inputDir + sysName + "/Load/ST.csv", ST_load);
	if (!status)	return false;
	
	status = readLoadData(inputDir + sysName + "/Load/RT.csv", RT_load);
	if (!status)	return false;
	
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
