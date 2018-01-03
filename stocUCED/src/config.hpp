//
//  config.h
//  Stoch-UC
//
//  Created by Semih Atakan on 12/8/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef config_h
#define config_h

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
	
};

const double pi = 3.1451;
const double loadShedPenaltyCoef = 1e5;
const double EPSzero = 1e-8;

const char delimiter = ',';

#endif /* config_h */
