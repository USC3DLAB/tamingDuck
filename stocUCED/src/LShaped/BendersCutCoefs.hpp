//
//  BendersCutCoefs.hpp
//  tamingDuck
//
//  Created by Semih Atakan on 2/17/18.
//  Copyright © 2018 Semih Atakan. All rights reserved.
//

#ifndef BendersCutCoefs_h
#define BendersCutCoefs_h

#include <vector>
#include "../misc.hpp"

/* Benders' Cuts */
struct BendersCutCoefs {
	double pi_b;
	vector< vector<double> > pi_T;
	
	void initialize(int numRows, int numCols) {
		resize_matrix(pi_T, numRows, numCols);
	};
	
	void reset() {
		pi_b = 0;
		for (int g=0; g < (int) pi_T.size(); g++) fill(pi_T[g].begin(), pi_T[g].end(), 0.0);
	};
};


#endif /* BendersCutCoefs_h */
