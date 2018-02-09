#include "instance.hpp"
#include "./stocProcess/stoc.hpp"

extern runType runParam;

instance::instance () {}

void instance::simulateScenarios(int numScen, bool fitModel) {
	int maxLookAhead = max(runParam.ED_numPeriods-1, runParam.ST_numPeriods-1);
	int numSimLengthInDays = ceil( (double)(runParam.numPeriods+maxLookAhead)/(60.0/runParam.baseTime) );
	
	simulations = createScenarioList(R, fitModel, powSys->path, stocElems, numSimLengthInDays, numScen);
	simulations.name = "ForecastData_simulations";
}

bool instance::initialize(PowSys *powSys, StocProcess *stoc, string RScriptsPath) {

	/* Initialize instance elements */
	this->powSys    = powSys;
	this->stoc	    = stoc;
	this->RScriptsPath = RScriptsPath;
	
	DAUC_t.resize( runParam.DA_numSolves );
	STUC_t.resize( runParam.DA_numSolves * runParam.ST_numSolves );
	ED_t.resize  ( runParam.DA_numSolves * runParam.ST_numSolves * runParam.ED_numSolves );
	
	detElems	= {"Load"};
	stocElems	= {"Solar", "Wind"};
	elems.resize( detElems.size() + stocElems.size() );
	merge( detElems.begin(), detElems.end(), stocElems.begin(), stocElems.end(), elems.begin() );
	
	/* Allocate memory for solution matrices */
	solution.allocateMem(powSys->numGen, runParam.numPeriods, powSys->numBus);

	/* Setup R environment */
	R["RScriptsPath"] = RScriptsPath;
	R["dataFolder"]   = powSys->path;	// set the data folder
	
	/* Chop off the Provided Time Series */
	int maxLookAhead = max(runParam.ED_numPeriods-1, runParam.ST_numPeriods-1);
	int numSimLengthInDays = ceil( (double)(runParam.numPeriods+maxLookAhead)/(60.0/runParam.baseTime) );

	vector<string> fileNames = {"DA", "RT"};
	for (int f=0; f<fileNames.size(); f++) {
		vector<int> stocIndices;
		for (auto it=stoc->mapTypeToIndex[fileNames[f]].begin(); it != stoc->mapTypeToIndex[fileNames[f]].end(); ++it) {
			if (find(elems.begin(), elems.end(), stoc->sp[*it].name) != elems.end()) {	// find the stoc process, and push its index to the list
				stocIndices.push_back(*it);
			}
		}
		observations.insert(pair<string, ScenarioType> (fileNames[f],
														createScenarioList(stoc, stocIndices, numSimLengthInDays, runParam.numPeriods/(60.0/runParam.baseTime), runParam.numRep)
														)
							);
		observations[fileNames[f]].name = fileNames[f] + "_observations";
	}
	
	/* simulate scenarios */
	simulateScenarios(1, true);			// no sampling, but only model fitting in here. 1 is set to avoid bugs.
	

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
