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
	int dataPeriodLengthInMins = runParam.baseTime;
	
	{/* actual observations */
		string filename = "RT";
		vector<int> idx = filename2idx(filename);
		
		// initialize the data structures
		actuals = createScenarioList(stoc, idx, numSimLengthInDays, runParam.numPeriods/(60.0/runParam.baseTime), runParam.numRep, dataPeriodLengthInMins);
		actuals.name = filename + "_obs";
		
		// post-processing
		boostRenewableSupply(actuals);
		correctSupplyExceedingCapacity(actuals);
	}
	
//	// print out summary
//	cout << actuals.name << ":" << endl;
//	int i=0;
//	for (auto it=actuals.mapVarNamesToIndex.begin(); it!=actuals.mapVarNamesToIndex.end(); ++it, ++i) {
//		cout << setw(4) << left << i << setw(30) << it->first << " (" << setprecision(2) << runParam.renewableCoef << "x)\t[OK]" << endl;
//	}
//	cout << endl;

	{	/* forecasts */
		vector<string> filenames = {"DA", "4HA"};

		for (int f=0; f<filenames.size(); f++) {
			string filename = filenames[f];
			vector<int> idx = filename2idx(filename);

			// initialize the data structures
			meanForecast.insert( pair<string, ScenarioType> (filename,
															 createScenarioList(stoc, idx, numSimLengthInDays, runParam.numPeriods/(60.0/runParam.baseTime), runParam.numRep, dataPeriodLengthInMins)) );
			meanForecast[filename].name = filename + "_forecasts";
			
			// post-processing
			boostRenewableSupply(meanForecast[filename]);
			correctSupplyExceedingCapacity(meanForecast[filename]);
			
			//	// print out summary
			//	cout << meanForecast[filename].name << ":" << endl;
			//	i=0;
			//	for (auto it=meanForecast[filename].mapVarNamesToIndex.begin(); it!=meanForecast[filename].mapVarNamesToIndex.end(); ++it, ++i) {
			//		cout << setw(4) << left << i << setw(30) << it->first << " (" << setprecision(2) << runParam.renewableCoef << "x)\t[OK]" << endl;
			//	}
			//	cout << endl;
		}
	}
	
//		The assumption here was a single forecast mechanism, that was updated during the day.
//		the "RT" forecast meant to be updated versions of the "DA" forecast.
//		We can keep with this tradition and only update "4HA" forecasts during the day, we
//		don't need "DA" forecasts during the day anyways.
//		// if updated forecasts are requested, initialize its structure
//		if (runParam.updateForecasts) {
//			filename = "RT";
//			meanForecast.insert( pair<string, ScenarioType> (filename, meanForecast["DA"]) );
//			meanForecast[filename].name = filename + "_forecasts";
//		}

	summary();

	return true;
}//END instance()

vector<int> instance::filename2idx(string filename) {
	vector<int> idx;
	
	// find the indices
	for (auto it=stoc->mapTypeToIndex[filename].begin(); it != stoc->mapTypeToIndex[filename].end(); ++it) {
		if (find(elems.begin(), elems.end(), stoc->sp[*it].name) != elems.end()) {	// find the stoc process, and push its index to the list
			idx.push_back(*it);
		}
	}
	
	return idx;
}

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
	// clear the previous info, if there is any
	simulations.clear();	// TODO: This clears the deterministic elements (load) as well.
	
	int maxLookAhead = max(runParam.ED_numPeriods-1, runParam.ST_numPeriods-1);
	int numSimLengthInDays = ceil( (double)(runParam.numPeriods+maxLookAhead)/(60.0/runParam.baseTime) );
	
#ifdef SAMPLE_USING_R
	simulations.insert( pair<string, ScenarioType> ("DA",
													createScenarioList(R, fitModel, powSys->path, stocElems, numSimLengthInDays, numScen, rep)) );
#else
	string simulationsFolder = "./Simulations/";
	simulations.insert( pair<string, ScenarioType> ("DA",
													createScenarioList(simulationsFolder, stocElems, numSimLengthInDays, numScen, rep)) );
#endif
	simulations["DA"].name = "DA_simulations";
	
	// post processing
	boostRenewableSupply(simulations["DA"]);
	correctSupplyExceedingCapacity(simulations["DA"]);
	
	// updated forecasts
	if (runParam.updateForecasts) {
		simulations.insert( pair<string, ScenarioType> ("RT", simulations["DA"]) );
		simulations["RT"].name = "RT_simulations";
	}
}

void instance::updateForecasts(int rep, int beginMin, int endMin) {
	
	// initializations
	const vector<double> weights = {
		0.900, 0.875, 0.850, 0.825,
		0.800, 0.775, 0.750, 0.725,
		0.700, 0.675, 0.650, 0.625,
		0.600, 0.575, 0.550, 0.525,
		0.500, 0.500, 0.500, 0.500};
	
	int DAIdx, RTIdx;
	int beg = beginMin/runParam.baseTime;	// current period
	int end = endMin/runParam.baseTime;		// period until updates will be done
	
	if (beg <= 0) return;	// no updates when no realization is available
//	if (end >= runParam.numPeriods) end = runParam.numPeriods-1;
	
	ScenarioType* DAForecasts = &meanForecast["DA"];
	ScenarioType* RTForecasts = &meanForecast["RT"];

	ScenarioType* DASimulations = simulations.size() > 0 ? &simulations["DA"] : NULL;
	ScenarioType* RTSimulations = simulations.size() > 0 ? &simulations["RT"] : NULL;
	
	// forecast updates
	for (auto it = actuals.mapVarNamesToIndex.begin(); it != actuals.mapVarNamesToIndex.end(); ++it) {
		// make the updates
		for (int t=beg; t<=end; t++) {			
			// get the indices associated with it->second (might be different if maps are filled in diff order)
			DAIdx = DAForecasts->mapVarNamesToIndex[ it->first ];
			RTIdx = RTForecasts->mapVarNamesToIndex[ it->first ];
			
			// point-forecasts
			RTForecasts->vals[rep][t][RTIdx] = actuals.vals[rep][beg-1][it->second] + DAForecasts->vals[rep][t][DAIdx] - DAForecasts->vals[rep][t-1][DAIdx];
			RTForecasts->vals[rep][t][RTIdx] = RTForecasts->vals[rep][t][RTIdx] * weights[t-beg] + DAForecasts->vals[rep][t][DAIdx] * (1-weights[t-beg]);
			
			if (RTForecasts->vals[rep][t][RTIdx] < 0)	RTForecasts->vals[rep][t][RTIdx] = 0.0;
			
			// simulations
			if (simulations.size() > 0) {
				// get the correct indices
				DAIdx = DASimulations->mapVarNamesToIndex[ it->first ];
				RTIdx = RTSimulations->mapVarNamesToIndex[ it->first ];
				
				// update the simulations
				for (int s=0; s<runParam.numTotScen; s++) {
					for (int t=beg; t<=end; t++) {
						RTSimulations->vals[s][t][RTIdx] = actuals.vals[rep][beg-1][it->second] + DASimulations->vals[s][t][DAIdx] - DASimulations->vals[s][t][DAIdx];
						RTSimulations->vals[rep][t][RTIdx] = RTSimulations->vals[rep][t][RTIdx] * weights[t-beg] + DASimulations->vals[rep][t][DAIdx] * (1-weights[t-beg]);
						
						if (RTSimulations->vals[rep][t][RTIdx] < 0)	RTSimulations->vals[rep][t][RTIdx] = 0.0;
					}
				}
			}
		}
	}
	
	
//	ScenarioType* DAForecasts = &meanForecast["DA"];
//	ScenarioType* RTForecasts = &meanForecast["RT"];
//
//	ScenarioType* DASimulations = simulations.size() > 0 ? &simulations["DA"] : NULL;
//	ScenarioType* RTSimulations = simulations.size() > 0 ? &simulations["RT"] : NULL;
//
//	int beg = beginMin/runParam.baseTime;	// current period
//	int end = endMin/runParam.baseTime;		// period until updates will be done
//
//	int DAIdx, RTIdx;
//	for (auto it = actuals.mapVarNamesToIndex.begin(); it != actuals.mapVarNamesToIndex.end(); ++it) {
//		// get the correct indices
//		DAIdx = DAForecasts->mapVarNamesToIndex[ it->first ];
//		RTIdx = RTForecasts->mapVarNamesToIndex[ it->first ];
//
//		// update point-forecasts
//		for (int t=beg; t<=end; t++) {
//			RTForecasts->vals[rep][t][RTIdx] = actuals.vals[rep][t][it->second] * futureInfoWeights[t-beg] + DAForecasts->vals[rep][t][DAIdx] * (1-futureInfoWeights[t-beg]);
//		}
//
//		// update simulations
//		if (simulations.size() > 0) {
//			// get the correct indices
//			DAIdx = DASimulations->mapVarNamesToIndex[ it->first ];
//			RTIdx = RTSimulations->mapVarNamesToIndex[ it->first ];
//
//			// update the simulations
//			for (int s=0; s<runParam.numTotScen; s++) {
//				for (int t=beg; t<=end; t++) {
//					RTSimulations->vals[s][t][RTIdx] = actuals.vals[rep][t][it->second] * futureInfoWeights[t-beg] + DASimulations->vals[s][t][DAIdx] * (1-futureInfoWeights[t-beg]);
//				}
//			}
//		}
//	}
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

void instance::boostRenewableSupply(ScenarioType &scenarioType) {
	// renewable-supply boost
	if ( fabs(runParam.renewableCoef-1) > EPSzero ) {
		for (auto it = powSys->generators.begin(); it != powSys->generators.end(); ++it) {
			if (it->type == Generator::SOLAR || it->type == Generator::WIND) {
				int genidx = scenarioType.mapVarNamesToIndex[it->name];	// get generator idx
				
				for (int r=0; r<scenarioType.vals.size(); r++) {			// increase the production in every scenario
					for (int t=0; t<scenarioType.vals[0].size(); t++) {
						scenarioType.vals[r][t][genidx] *= runParam.renewableCoef;
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

	/* summarized statistics */
	vector<string> fields;
	fields.push_back("#ofOnGen");
	fields.push_back("UsedGen");
	fields.push_back("OverGen");
	fields.push_back("UsedGenCost");
	fields.push_back("OverGenCost");
	fields.push_back("NoLoadCost");
	fields.push_back("StartUpCost");
	fields.push_back("LoadShed");
	fields.push_back("CO2(kg)");
	fields.push_back("NOX(kg)");
	fields.push_back("SO2(kg)");
	fields.push_back("Capacity");
	fields.push_back("RampUpCap");
	fields.push_back("RampDownCap");
	fields.push_back("MinGenReq");

	/* print solutions */
	ofstream output;
	status = open_file(output, filepath + "_commitments.sol");
	if (!status) goto finalize;

	print_matrix(output, solution.x, delimiter, 0);
	output.close();

	status = open_file(output, filepath + "_genSTUC.sol");
	if (!status) goto finalize;

	print_matrix(output, solution.g_STUC, delimiter, 2);
	output.close();

	status = open_file(output, filepath + "_genDAUC.sol");
	if (!status) goto finalize;
	
	print_matrix(output, solution.g_DAUC, delimiter, 2);
	output.close();

	status = open_file(output, filepath + "_genED.sol");
	if (!status) goto finalize;

	print_matrix(output, solution.g_ED, delimiter, 2);
	output.close();

	status = open_file(output, filepath + "_loadShedED.sol");
	if (!status) goto finalize;
	
	print_matrix(output, solution.loadShed_ED, delimiter, 2);
	output.close();
	
	status = open_file(output, filepath + "_usedGenED.sol");
	if (!status) goto finalize;
	
	print_matrix(output, solution.usedGen_ED, delimiter, 2);
	output.close();

	status = open_file(output, filepath + "_overGenED.sol");
	if (!status) goto finalize;
	
	print_matrix(output, solution.overGen_ED, delimiter, 2);
	output.close();

	/* print summarized statistics */
	status = open_file(output, filepath + "_stats.sol");
	if (!status) goto finalize;
	
	// headers
	output << "Time";
	for (int f=0; f<fields.size(); f++) {
		output << "\t" << fields[f];
	}
	output << endl;

	if(runParam.useGenHistory) {
		cout << "Start up costs may be underestimated" << endl;
	}

	// data
	for (int t=0; t<runParam.numPeriods; t++) {
		output << setfill('0') << setw(2) << timeInfo->tm_hour << ":" << setfill('0') << setw(2) << timeInfo->tm_min << "\t";
		
		double coef = runParam.ED_resolution/60.0;
		
		// reset stats to 0.0
		map<string, double> stats;
		for (int f=0; f<fields.size(); f++) {
			stats.insert(pair<string,double> (fields[f], 0.0));
		}
		
		for (int g=0; g<powSys->numGen; g++) {
			// open generator count
			stats["#ofOnGen"] += solution.x[g][t];
			
			// capacities
			stats["Capacity"] += powSys->generators[g].maxCapacity * solution.x[g][t];
			stats["RampUpCap"] += powSys->generators[g].rampUpLim * solution.x[g][t];
			stats["RampDownCap"] += powSys->generators[g].rampDownLim * solution.x[g][t];
			stats["MinGenReq"] += powSys->generators[g].minGenerationReq * solution.x[g][t];
			
			// production amounts
			stats["UsedGen"] += solution.usedGen_ED[g][t];
			stats["OverGen"] += solution.overGen_ED[g][t];
			
			// costs
			stats["UsedGenCost"] += powSys->generators[g].variableCost * coef * solution.usedGen_ED[g][t];
			stats["OverGenCost"] += powSys->generators[g].variableCost * coef * solution.overGen_ED[g][t];
			stats["NoLoadCost"] += powSys->generators[g].noLoadCost   * coef * solution.x[g][t];
			
			if (t > 0) {
				stats["StartUpCost"] += powSys->generators[g].startupCost * max(solution.x[g][t]-solution.x[g][t-1], 0.0);
			} else if(runParam.useGenHistory) {
//				cout << "Start up costs may be underestimated" << endl;
			}
			
			// emissions
			stats["CO2(kg)"] += powSys->generators[g].CO2_emission_base * coef * solution.x[g][t];
			stats["CO2(kg)"] += powSys->generators[g].CO2_emission_var * coef * solution.g_ED[g][t];
			stats["NOX(kg)"] += powSys->generators[g].NOX_emission_base * coef * solution.x[g][t];
			stats["NOX(kg)"] += powSys->generators[g].NOX_emission_var * coef * solution.g_ED[g][t];
			stats["SO2(kg)"] += powSys->generators[g].SO2_emission_base * coef * solution.x[g][t];
			stats["SO2(kg)"] += powSys->generators[g].SO2_emission_var * coef * solution.g_ED[g][t];
		}
		
		// shed demand
		for (int b = 0; b < powSys->numBus; b++) {
			stats["LoadShed"] += solution.loadShed_ED[b][t];
		}
		
		// print
		for (int f=0; f<fields.size(); f++) {
			output << stats[fields[f]];
			if (f != fields.size()-1) output << "\t";
		}
		output << endl;
		
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

/****************************************************************************
 * addCurrentSolToSolList
 ****************************************************************************/
void instance::addCurrentSolToSolList() {
	solList.push_back(solution);
	out() << "Decisions are saved." << endl;
}

