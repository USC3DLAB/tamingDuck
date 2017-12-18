//
//  Generator.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef Generator_hpp
#define Generator_hpp

#include <stdio.h>
#include <string>

#include "Bus.hpp"

using namespace std;

class Bus;

class Generator {
    
public:
    Generator ();
    
    void setType (string typeName); // sets generator type given its name
    
    // generator identifiers
	string name;	// provided by the data
	int id;			// assigned by us
    
    // generator characteristics
    double maxCapacity;         // MW
    double minGenerationReq;    // MW
    
    double rampUpLim;           // MW / minute
    double rampDownLim;         // MW / minute
    
    double variableCost;        // $ / MWh
    double startupCost;         // $ / startup
    double noLoadCost;          // $ / hour
    
    int minUpTime;              // hours
    int minDownTime;            // hours
    
    bool isMustRun;             // 1 if generator must always run, and its production must always be consumed by the grid.
    bool isBaseLoadGen;         // 1 if generator must be scheduled in the long-term UC, 0 if it must be scheduled in short-term UC.
    
    string  connectedBusName;   // name of the bus that this generator is connected to
    Bus     *connectedBus;      // bus that this generator is connected to

    enum GeneratorType {
        BIOMASS,
        NATURALGAS,
        OIL,
        GEOTHERMAL,
        HYDRO,
        COAL,
        SOLAR,
        WIND,
        OTHER
    };
    GeneratorType type;
};

#endif /* Generator_hpp */
