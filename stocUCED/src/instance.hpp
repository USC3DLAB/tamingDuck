
/*************************************************
 
 Instance for a power production planning problem
 
 ************************************************/

#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include "config.hpp"
#include "./stocProcess/stoc.hpp"
#include "./powerSys/PowSys.hpp"
#include "solution.hpp"
#include "misc.hpp"

#include <list>

using namespace std;

class instance {
	
public:
    instance ();
	bool initialize			(PowSys *powSys, StocProcess *stoc, string RScriptsPath);
	bool printSolution		(string filepath);
	void simulateScenarios	(int numScen, bool fitModel, int rep);
	void setRSeed			(int rep);
	
	PowSys		*powSys;
	StocProcess	*stoc;
	Solution	solution;
	
	list<Solution> solList;
	void addCurrentSolToSolList ();
	
	vector<string> elems;			// list of all processes
	vector<string> detElems;		// list of deterministic elements in the model
	vector<string> stocElems;		// list of stochastic elements in the model
	
	map<string, ScenarioType> simulations;	// DA simulations, RT simulations (updated according to past actual-realizations)
	map<string, ScenarioType> meanForecast;	// DA mean-forecast, RT mean-forecast (updated according to past actual-realizations)
	ScenarioType actuals;					// actual realizations
	
	void updateForecasts(int rep, int beginMin, int endMin);
	
	/* log keeping */
	ofstream&	out ();
	string		getLogStreamName ();
	bool		openLogFile (string filename);
	void		closeLogFile();
	
private:
	void summary();
	
	/* log keeping */
	ofstream	log_stream;
	string		log_stream_name;
	
#ifdef SAMPLE_USING_R
	RInside R;
#endif
	string RScriptsPath;
	
	void correctSupplyExceedingCapacity(ScenarioType &scenarioType);
	void boostRenewableSupply(ScenarioType &scenarioType);
<<<<<<< HEAD
=======
	
	vector<int> filename2idx (string filename);
	
>>>>>>> 2247e8a849f9d0cf0fc44445ea459889ee1f793e
};

#endif
