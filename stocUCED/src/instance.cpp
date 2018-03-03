#include "instance.hpp"
#include "./stocProcess/stoc.hpp"

extern runType runParam;

instance::instance () {}

bool instance::initialize(PowSys *powSys, StocProcess *stoc, string RScriptsPath) {

	/* Initialize instance elements */
	this->powSys    = powSys;
	this->stoc	    = stoc;
	this->RScriptsPath = RScriptsPath;
	
	detElems	= {"Load"};
	stocElems	= {"Solar", "Wind"};
	elems.resize( detElems.size() + stocElems.size() );
	merge( detElems.begin(), detElems.end(), stocElems.begin(), stocElems.end(), elems.begin() );
	
	/* Allocate memory for solution matrices */
	solution.allocateMem(powSys->numGen, runParam.numPeriods, powSys->numBus);

#ifdef SAMPLE_USING_R
	/* Setup R environment */
	R["RScriptsPath"] = RScriptsPath;
	R["dataFolder"]   = powSys->path;	// set the data folder
#endif
	
	/* Chop off the Provided Time Series */
	int maxLookAhead = max(runParam.ED_numPeriods-1, runParam.ST_numPeriods-1);
	int numSimLengthInDays = ceil( (double)(runParam.numPeriods+maxLookAhead)/(60.0/runParam.baseTime) );

	vector<string> fileNames = {"DA", "RT"};
	for (int f = 0; f < (int) fileNames.size(); f++) {
		vector<int> stocIndices;
		for (auto it=stoc->mapTypeToIndex[fileNames[f]].begin(); it != stoc->mapTypeToIndex[fileNames[f]].end(); ++it) {
			if (find(elems.begin(), elems.end(), stoc->sp[*it].name) != elems.end()) {	// find the stoc process, and push its index to the list
				stocIndices.push_back(*it);
			}
		}
		
		//TODO: Generalize
		//int dataPeriodLengthInMins = fileNames[f] == "DA" ? 60 : 15;
		int dataPeriodLengthInMins = 15;
		
		observations.insert(pair<string, ScenarioType> (fileNames[f],
														createScenarioList(stoc, stocIndices, numSimLengthInDays, runParam.numPeriods/(60.0/runParam.baseTime), runParam.numRep, dataPeriodLengthInMins)));
		observations[fileNames[f]].name = fileNames[f] + "_observations";
		
		cout << observations[fileNames[f]].name << ":" << endl;
		int i=0;
		for (auto it=observations[fileNames[f]].mapVarNamesToIndex.begin(); it!=observations[fileNames[f]].mapVarNamesToIndex.end(); ++it) {
			cout << setw(4) << left << ++i << setw(30) << it->first << "\t[OK]" << endl;
		}
		
		correctSupplyExceedingCapacity(observations[fileNames[f]]);
	}
	
	/* simulate scenarios */
	simulateScenarios(1, true, 0);			// no sampling, but only model fitting in here. 1 is set to avoid bugs.
	

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

/****************************************************************************
 * simulateScenarios
 * Uses R scripts to simulate multiple time series
 ****************************************************************************/
void instance::simulateScenarios(int numScen, bool fitModel, int rep) {
	int maxLookAhead = max(runParam.ED_numPeriods-1, runParam.ST_numPeriods-1);
	int numSimLengthInDays = ceil( (double)(runParam.numPeriods+maxLookAhead)/(60.0/runParam.baseTime) );
	
#ifdef SAMPLE_USING_R
	simulations = createScenarioList(R, fitModel, powSys->path, stocElems, numSimLengthInDays, numScen, rep);
#else
	string simulationsFolder = "../Release/Simulations/";
	simulations = createScenarioList(simulationsFolder, stocElems, numSimLengthInDays, numScen, rep);
#endif
	
	
	simulations.name = "ForecastData_simulations";
	
	correctSupplyExceedingCapacity(simulations);
}

/****************************************************************************
 * correctSupplyExceedingCapacity
 * Both the simulations and the provided data may lead to supply amounts,
 * which exceed generator capacities. This function goes over the simulated
 * or the original time series in order to set all supply = min(supply, cap).
 ****************************************************************************/
void instance::correctSupplyExceedingCapacity(ScenarioType &scenarioType) {
	// iterate over generators
	for (auto genItr = powSys->generators.begin(); genItr != powSys->generators.end(); ++genItr) {
		auto it = scenarioType.mapVarNamesToIndex.find( (*genItr).name );
		
		// find a generator that is registered in the "scenarioType"
		if ( it != scenarioType.mapVarNamesToIndex.end() ) {
			for (int s=0; s< (int) scenarioType.vals.size(); s++) {
				for (int t=0; t< (int) scenarioType.vals[s].size(); t++) {
					
					// for all scenarios and time periods, make sure the capacity is not exceeded
					if ( scenarioType.vals[s][t][it->second] > (*genItr).maxCapacity ) {
						scenarioType.vals[s][t][it->second] = (*genItr).maxCapacity;
					}
				}
			}
		}
	}
}

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

	status = open_file(output, filepath + "_genDAUC.sol");
	if (!status) goto finalize;

	print_matrix(output, solution.g_DAUC, delimiter, 2);
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
		
		double coef = runParam.ED_resolution/60.0;
		
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
		for (int b = 0; b < powSys->numBus; b++) {
			stats[7] += solution.loadShed_ED[b][t];
		}
		
		for (int i = 0; i < (int) stats.size()-1; i++) {
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

void instance::setRSeed (int rep) {
#ifdef SAMPLE_USING_R
	R["rep"] = rep;
	R.parseEval("set.seed(rep)");
#endif
}

/****************************************************************************
 * out
 * returns a pointer to the log stream
 ****************************************************************************/
ofstream& instance::out() {
	return log_stream;
}

/****************************************************************************
 * getLogStreamName
 * returns the name of the log stream
 ****************************************************************************/
string instance::getLogStreamName() {
	return log_stream_name;
}

/****************************************************************************
 * openLogFile
 * opens a log file with the input name. 
 * log files can be used to print out optimization logs from CPLEX or SD.
 ****************************************************************************/
bool instance::openLogFile(string filename) {
	bool status = open_file(log_stream, filename);
	if (status) {
		log_stream_name = filename;
	}
	
	return status;
}

/****************************************************************************
 * closeLogFile
 ****************************************************************************/
void instance::closeLogFile() {
	log_stream.close();
}
