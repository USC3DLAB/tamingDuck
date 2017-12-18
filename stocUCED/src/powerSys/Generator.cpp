//
//  Generator.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "Generator.hpp"

Generator::Generator() {}

void Generator::setType(string typeName) {
    if (typeName.compare("Biomass") == 0)               type = BIOMASS;
    else if (typeName.compare("Natural gas") == 0)      type = NATURALGAS;
    else if (typeName.compare("Oil") == 0)              type = OIL;
    else if (typeName.compare("Geothermal") == 0)       type = GEOTHERMAL;
    else if (typeName.compare("Hydro") == 0)            type = HYDRO;
    else if (typeName.compare("Solar") == 0)            type = SOLAR;
    else if (typeName.compare("Coal") == 0)             type = COAL;
    else if (typeName.compare("Wind") == 0)             type = WIND;
    else if (typeName.compare("Other") == 0)            type = OTHER;
    else {
        perror( ("Error: Unknown generator type: " + typeName + "\n").c_str() );
    }
}
