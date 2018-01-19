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
	// SMH: How about the constructor? or just hardcode it in the header file? Seems like the 3-layer framework quite hardcoded anyways. 
	this->hierarchy = {"DA", "ST", "RT"};

	/* Allocate memory for solution matrices */
	solution.allocateMem(powSys->numGen, runParam.numPeriods);

	/* Setup all the deterministic observations */
	for ( int l = 0; l < (int) this->hierarchy.size(); l++ ) {
		
		vector<int> indices;
		for (auto it=stoc->mapTypeToIndex[this->hierarchy[l]].begin(); it != stoc->mapTypeToIndex[this->hierarchy[l]].end(); ++it) {
			if ( find(detElems.begin(), detElems.end(), stoc->sp[*it].name) != detElems.end() ) {	// if the name is found, push it into the list of indices
				indices.push_back(*it);
			}
		}
		// TODO: this num rep issue will cause trouble
		ScenarioType temp = createScenarioList(stoc, indices, runParam.numPeriods+runParam.ED_numPeriods-1, runParam.numPeriods, runParam.numRep);
		temp.name = this->hierarchy[l] + "_DetObs";
		this->detObserv.push_back(temp);
	}

	/* Setup all the stochastic observations */
	for ( int l = 0; l < (int) this->hierarchy.size(); l++ ) {
		
		vector<int> indices;
		for (auto it=stoc->mapTypeToIndex[this->hierarchy[l]].begin(); it != stoc->mapTypeToIndex[this->hierarchy[l]].end(); ++it) {
			cout << stoc->sp[*it].name << " " << this->hierarchy[l] << " " << *it << endl;
			if ( find(stocElems.begin(), stocElems.end(), stoc->sp[*it].name) != stocElems.end() ) {	// if the name is found, push it into the list of indices
				indices.push_back(*it);
			}
		}
		// TODO: this num rep issue will cause trouble
		ScenarioType temp = createScenarioList(stoc, indices, runParam.numPeriods+runParam.ED_numPeriods-1, runParam.numPeriods, runParam.numRep);
		temp.name = this->hierarchy[l] + "_StochObs";
		this->stocObserv.push_back(temp);
	}

	summary();

	return true;
}//END instance()


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
