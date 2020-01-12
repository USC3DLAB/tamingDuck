//
//  solution.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 11/11/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef solution_hpp
#define solution_hpp

#include <stdio.h>
#include <vector>

#include "misc.hpp"

struct Solution {
    
    Solution () {}
    
    void allocateMem (int numGen, int periods, int numBus, int numBatteries) {
        resize_matrix(x,	numGen, periods);
		resize_matrix(g_DAUC, numGen, periods);
		resize_matrix(g_STUC, numGen, periods);
        resize_matrix(g_ED, numGen, periods);
		resize_matrix(overGen_ED, numGen, periods);
		resize_matrix(usedGen_ED, numGen, periods);
		resize_matrix(loadShed_ED, numBus, periods);
		resize_matrix(btState_ED, numBatteries, periods);
		resize_matrix(btState_UC, numBatteries, periods);
		resize_matrix(btFlow_ED, numBatteries, periods);
    }
    
    vector< vector<double> > x, g_UC, g_ED, overGen_ED, usedGen_ED, loadShed_ED, btState_ED, btState_UC, btFlow_ED;
};

#endif /* solution_hpp */
