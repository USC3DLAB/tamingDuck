//
//  Line.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef Line_hpp
#define Line_hpp

#include <stdio.h>
#include <string>

#include "Bus.hpp"

class Bus;

class Line {

public:
    Line ();
    
    // line identifiers
    string name;
    
    // line characteristics
    double minFlowLim;  // MWh
    double maxFlowLim;  // MWh
    double susceptance; // ?
    
    Bus *orig;  // origin bus
    Bus *dest;  // destination bus
    
};
#endif /* Line_hpp */
