#include "instance.hpp"
#include "misc.hpp"

// TODO: SMH: does this need to be in here?
extern runType runParam;

instance::instance () {}

void instance::initialize(PowSys *powSys, scenarios *stoc, string inputDir) {
	this->powSys = powSys;
	this->stoc	 = stoc;
	
	solution.allocateMem(powSys->numGen, (int)round(runParam.horizon/runParam.ED_resolution));

	readLoadData(inputDir + "Load/DA.csv", DA_load);
	readLoadData(inputDir + "Load/ST.csv", ST_load);
	readLoadData(inputDir + "Load/RT.csv", RT_load);
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
