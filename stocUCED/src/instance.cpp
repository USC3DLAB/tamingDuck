#include "instance.hpp"
#include "./stocProcess/stoc.hpp"

extern runType runParam;

instance::instance () {}

bool instance::initialize(PowSys *powSys, StocProcess *stoc, vector<string> stocElems, vector<string> detElems) {

	/* Initialize instance elements */
	this->powSys    = powSys;
	this->stoc	    = stoc;
	this->detElems  = detElems;
	this->stocElems = stocElems;

	// TODO: Move this initialization to a better place.
	// SMH: How about the constructor? or just hardcode it in the header file?
	this->hierarchy = {"DA", "ST", "RT"};

	/* Allocate memory for solution matrices */
	solution.allocateMem(powSys->numGen, runParam.numPeriods);

	/* Setup all the deterministic observations */
	for ( int l = 0; l < (int) this->hierarchy.size(); l++ ) {
		for ( int n = 0; n < (int) detElems.size(); n++ ) {
			vector<int> indices;
			for (auto it=stoc->mapTypeToIndex[this->hierarchy[l]].begin(); it != stoc->mapTypeToIndex[this->hierarchy[l]].end(); ++it) {
				if ( find(detElems.begin(), detElems.end(), stoc->sp[*it].name) != detElems.end() ) {	// if the name is found, push it into the list of indices
					indices.push_back(*it);
				}
			}
			ScenarioType temp = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);
			this->detObserv.push_back(temp);
		}
	}

	/* Setup all the stochastic observations */
	for ( int l = 0; l < (int) this->hierarchy.size(); l++ ) {
		for ( int n = 0; n < (int) stocElems.size(); n++ ) {
			vector<int> indices;
			for (auto it=stoc->mapTypeToIndex[this->hierarchy[l]].begin(); it != stoc->mapTypeToIndex[this->hierarchy[l]].end(); ++it) {
				if ( find(stocElems.begin(), stocElems.end(), stoc->sp[*it].name) != stocElems.end() ) {	// if the name is found, push it into the list of indices
					indices.push_back(*it);
				}
			}
			ScenarioType temp = createScenarioList(stoc, indices, runParam.numPeriods, runParam.numRep);
			this->stocObserv.push_back(temp);
		}
	}

	summary();

	return true;
}//END instance()

/*
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
 */

void instance::summary() {
	cout << "------------------------------------------------------------------" << endl;
	printf("%-23s%s%s\n", "Power System", ": ", powSys->name.c_str());
	if ( !detElems.empty() ) {
		printf("%-23s%s", "Deterministic elements", ": ");
		for (auto i = detElems.begin(); i != detElems.end(); ++i)
			std::cout << *i << ' ';
		cout << endl;
	}
	else {
 		printf("%-23s%s%s\n", "Deterministic elements", ": ", "None");
	}
	
	if ( !stocElems.empty() ) {
		printf("%-23s%s", "Stochastic elements", ": ");
		for (auto i = stocElems.begin(); i != stocElems.end(); ++i)
			std::cout << *i << ' ';
		cout << endl;
	}
	else {
		printf("%-23s%s%s\n", "Stochastic elements", ": ", "None");
	}
	cout << "------------------------------------------------------------------" << endl;

}// summary()

bool instance::printSolution(string filepath) {
	bool status;
	
	ofstream output;
	status = open_file(output, filepath + "_commitments.sol");
	if (!status) goto finalize;
	
	print_matrix(output, solution.x, delimiter, 0);
	output.close();
	
	status = open_file(output, filepath + "_genUC.sol");
	if (!status) goto finalize;
	
	print_matrix(output, solution.g_UC, delimiter, 2);
	output.close();
	
	status = open_file(output, filepath + "_genED.sol");
	if (!status) goto finalize;
	
	print_matrix(output, solution.g_ED, delimiter, 2);
	output.close();
	
finalize:
	return status;
}
