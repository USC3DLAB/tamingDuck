//
//  Battery.hpp
//  tamingDuck
//
//  Created by Semih Atakan on 1/20/19.
//  Copyright © 2019 Semih Atakan. All rights reserved.
//

#ifndef Battery_hpp
#define Battery_hpp

#include <stdio.h>
#include "Bus.hpp"
using namespace std;

class Bus;
class Battery {
public:
	Battery ();

	// battery identifiers
	string name;	// provided by the data
	int id;			// assigned by us
	
	// battery characteristics
	double maxCapacity;         // MW
	double dissipationCoef;		// dissipation coef per hr (if capacity is I_t at t, 1 hr later, it will be I_t * dissipationCoef)
	double conversionLossCoef;	// conversion loss per hr due to charging / discharging
	string  connectedBusName;   // name of the bus that this generator is connected to
	Bus     *connectedBus;      // bus that this generator is connected to
};
#endif /* Battery_hpp */
