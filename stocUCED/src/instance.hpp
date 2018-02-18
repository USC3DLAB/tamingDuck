
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
	
	vector<string> elems;			// list of all processes
	vector<string> detElems;		// list of deterministic elements in the model
	vector<string> stocElems;		// list of stochastic elements in the model
	
	ScenarioType simulations;
	
	map< string, ScenarioType > observations;
	
	
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
};

#endif
