//
//  config.h
//  Stoch-UC
//
//  Created by Semih Atakan on 12/8/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef config_h
#define config_h

struct runType {
	double	horizon;		// Total horizon in minutes
	double	DA_horizon;		// in minutes
	double	ST_horizon;		// in minutes
	double	ED_horizon;		// in minutes
	int		DA_numPeriods;	// Number of periods for a given frequency
	int		ST_numPeriods;
	int		ED_numPeriods;
	double	DA_resolution;	// in minutes
	double	ST_resolution;	// in minutes
	double	ED_resolution;	// in minutes
	double	DA_frequency;	// in minutes
	double	ST_frequency;	// in minutes
	double	ED_frequency;	// in minutes
	int 	DA_numSolves;	// Number of times DA problem is solved
	int		ST_numSolves;	// Number of times ST problem needs to be solved within the horizon of DA
	int		ED_numSolves;	// Number of times ED problem needs to be solved within the horizon of ED
};

enum ProblemType {
	DayAhead,
	ShortTerm
};

enum ModelType {
	System,
	Transmission
};

const double pi = 3.1451;
const double loadShedPenaltyCoef = 1e5;

const char delimiter = ',';

#endif /* config_h */
