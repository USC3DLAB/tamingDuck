//
//  Bus.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef Bus_hpp
#define Bus_hpp

#include <stdio.h>
#include <string>
#include <vector>

#include "Generator.hpp"

using namespace std;

class Generator;

class Bus {

public:
    Bus ();

    // bus identifiers
    int id;
    string name;

    const double pi = 3.1451;
    // TODO: duplicate, another copy is in intance.hpp
    
    // bus characteristics
    int     regionId;

    double  loadPercentage;
    double  maxPhaseAngle =  pi;
    double  minPhaseAngle = -pi;
    
    vector<Generator*> connectedGenerators;
};

#endif /* Bus_hpp */
