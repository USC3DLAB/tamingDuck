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
	solution.allocateMem(powSys->numGen, runParam.numPeriods, powSys->numBus);

	int lookahead = max(runParam.ED_numPeriods-1, runParam.ST_numPeriods-1);
	/* Setup all the deterministic observations */
	for ( int l = 0; l < (int) this->hierarchy.size(); l++ ) {
		
		vector<int> indices;
		for (auto it=stoc->mapTypeToIndex[this->hierarchy[l]].begin(); it != stoc->mapTypeToIndex[this->hierarchy[l]].end(); ++it) {
			if ( find(detElems.begin(), detElems.end(), stoc->sp[*it].name) != detElems.end() ) {	// if the name is found, push it into the list of indices
				indices.push_back(*it);
			}
		}
		ScenarioType temp = createScenarioList(stoc, indices, runParam.numPeriods+lookahead, runParam.numPeriods, runParam.numRep);
		temp.name = this->hierarchy[l] + "_DetObs";
		this->detObserv.push_back(temp);
	}

	/* Setup all the stochastic observations */
	for ( int l = 0; l < (int) this->hierarchy.size(); l++ ) {
		
		vector<int> indices;
		for (auto it=stoc->mapTypeToIndex[this->hierarchy[l]].begin(); it != stoc->mapTypeToIndex[this->hierarchy[l]].end(); ++it) {
			if ( find(stocElems.begin(), stocElems.end(), stoc->sp[*it].name) != stocElems.end() ) {	// if the name is found, push it into the list of indices
				indices.push_back(*it);
			}
		}
		/*******************************************************************************
		 Important: The created scenario list will include more fields than the number
		 of periods, so that the ED problem can be solved in a rolling horizon fashion. 
		 Accordingly, for instance, for a 24-hour problem, a set of realizations must 
		 include 24*ED_resolution + runParam.ED_numPeriods-1 elements
		 *******************************************************************************/
		ScenarioType temp = createScenarioList(stoc, indices, runParam.numPeriods+lookahead, runParam.numPeriods, runParam.numRep);
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
	
	/* time related */
	time_t rawTime;
	struct tm * timeInfo;
	time ( &rawTime );
	timeInfo = localtime( &rawTime );
	timeInfo->tm_min = 0;	timeInfo->tm_hour = 0; mktime(timeInfo);

	/* print solutions */
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

	status = open_file(output, filepath + "_stats.sol");
	if (!status) goto finalize;
	output << "Time\t#ofOnGen\tUsedGen\tOverGen\tUsedGenCost\tOverGenCost\tNoLoadCost\tStartUpCost\tLoadShed\tCO2(kg)\tNOX(kg)\tSO2(kg)" << endl;
	for (int t=0; t<runParam.numPeriods; t++) {
		output << setfill('0') << setw(2) << timeInfo->tm_hour << ":" << setfill('0') << setw(2) << timeInfo->tm_min << "\t";
		
		double coef = 60.0/runParam.ED_resolution;
		
		vector<double> stats (11, 0.0);
		for (int g=0; g<powSys->numGen; g++) {
			// open generator count
			stats[0] += solution.x[g][t];
			
			// production amounts
			stats[1] += solution.usedGen_ED[g][t];
			stats[2] += solution.overGen_ED[g][t];
			
			// costs
			stats[3] += powSys->generators[g].variableCost * coef * solution.usedGen_ED[g][t];
			stats[4] += powSys->generators[g].variableCost * coef * solution.overGen_ED[g][t];
			stats[5] += powSys->generators[g].noLoadCost   * coef * solution.x[g][t];
			
			if (t > 0) {
				stats[6] += powSys->generators[g].startupCost * max(solution.x[g][t]-solution.x[g][t-1], 0.0);
			}
			
			// emissions
			stats[8] += powSys->generators[g].CO2_emission_base * coef * solution.x[g][t];
			stats[8] += powSys->generators[g].CO2_emission_var * coef * solution.g_ED[g][t];
			stats[9] += powSys->generators[g].NOX_emission_base * coef * solution.x[g][t];
			stats[9] += powSys->generators[g].NOX_emission_var * coef * solution.g_ED[g][t];
			stats[10] += powSys->generators[g].SO2_emission_base * coef * solution.x[g][t];
			stats[10] += powSys->generators[g].SO2_emission_var * coef * solution.g_ED[g][t];
		}
		
		// shed demand
		for (int b=0; b<powSys->numBus; b++) {
			stats[7] += solution.loadShed_ED[b][t];
		}
		
		for (int i=0; i<stats.size()-1; i++) {
			output << stats[i] << "\t";
		}
		output << stats[ stats.size()-1 ] << endl;
		
		timeInfo->tm_min += runParam.baseTime;
		mktime(timeInfo);
	}
	
	output.close();
	
	finalize:
	return status;
}
