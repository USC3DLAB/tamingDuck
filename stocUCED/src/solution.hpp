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

/* State variable which provides the status of generator.
 *
 * s      = start up
 * x      = state of the generator, 1 if operational and 0 otherwise
 * z      = shut down
 * gInit  = initial generation
 *
 */
struct Solution {
    
    Solution () {}
    
    void allocateMem ( int numGen, int periods ) {
        resize_matrix(x, numGen, periods);
        resize_matrix(g_UC, numGen, periods);
        resize_matrix(g_ED, numGen, periods);
    }
    
    vector< vector<double> > x, g_UC, g_ED;
};

#endif /* solution_hpp */
