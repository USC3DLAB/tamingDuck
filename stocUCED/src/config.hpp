//
//  config.h
//  Stoch-UC
//
//  Created by Semih Atakan on 12/8/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef config_h
#define config_h

// MACROS
// #define SAMPLE_USING_R
#define BOOST_PARALLEL_LIBS		// Parallel programming libraries of boost is being used

enum ProblemType {
	DayAhead,
	ShortTerm
};

enum ModelType {
	System,
	Transmission
};

enum LevelType {
	DA,
	ST,
	RT
};

enum SettingType {
	DUC_DED,
	DUC_SED,
	SUC_SED
};

enum Setting {
	DETERMINISTIC,
	STOCHASTIC,
	NA
};

struct runType {
	double 	baseTime;		// in minutes, the highest resolution period length
	double	horizon;		// Total horizon in minutes
	double	DA_horizon;		// in minutes
	double	ST_horizon;		// in minutes
	double	ED_horizon;		// in minutes
	int		DA_numPeriods;	// Number of periods for a given frequency
	int		ST_numPeriods;
	int		ED_numPeriods;
	int		numPeriods;		// Total # of periods
	double	DA_resolution;	// in minutes
	double	ST_resolution;	// in minutes
	double	ED_resolution;	// in minutes
	double	DA_frequency;	// how frequently the DA problem is solved (in minutes)
	double	ST_frequency;	// how frequently the ST problem is solved (in minutes)
	double	ED_frequency;	// how frequently the ED problem is solved (in minutes)
	int 	DA_numSolves;	// Number of times DA problem is solved
	int		ST_numSolves;	// Number of times ST problem needs to be solved within the horizon of DA
	int		ED_numSolves;	// Number of times ED problem needs to be solved within the horizon of ED
	
	SettingType settingType;
	
	int		numRep;			// Number of replications
	int		numTotScen;		// Number of time-series generated in total
	int		numLSScen;		// Number of time-series needed to sample for L-Shaped
	
	int 	startRep;		// which replication to start with (minimum is 1)
	
	bool 	useGenHistory;			// true if generator histories from the previous run is being used, false oth.
	double 	spinResPerc;			// spinning reserve percentage
	double	renewableMultiplier;	// renewable-supply will be multiplied with this factor
	bool	updateForecasts;		// real-time forecast updates
};

const double pi = 3.1451;
const double loadShedPenaltyCoef = 5000;
const double overGenPenaltyCoef = 25.0;
const double EPSzero = 1e-8;

const char delimiter = ',';

// Tolerances
const double AbsLShapedCutOptTol = 1;
const double RelLShapedCutOptTol = 1e-5;

// Thread management
const unsigned short LShapedMasterCPXThreads = 1;
const unsigned short LShapedSubprobCPXThreads = 1;
#ifdef BOOST_PARALLEL_LIBS
#include <thread>
const unsigned short LShapedSubprobThreads = std::thread::hardware_concurrency(); //boost::thread::hardware_concurrency();
#else
const unsigned short LShapedSubprobThreads = 1;
#endif

#endif /* config_h */
