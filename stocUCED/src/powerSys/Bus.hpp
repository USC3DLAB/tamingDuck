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
#include "Battery.hpp"
#include "../config.hpp"

using namespace std;

class Generator;
class Battery;

class Bus {

public:
    Bus ();

    // bus identifiers
	string name;	// provided by the data
	int id;			// assigned by us

    // bus characteristics
    int     regionId;

    double  loadPercentage;
    double  maxPhaseAngle =  pi;
    double  minPhaseAngle = -pi;
    
    vector<Generator*> 	connectedGenerators;
	vector<Battery*>	connectedBatteries;
};

#endif /* Bus_hpp */
