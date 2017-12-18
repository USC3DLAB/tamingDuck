//
//  PowSys.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef PowSys_hpp
#define PowSys_hpp

#include <stdio.h>
#include <vector>
#include <fstream>
#include <string>
#include <map>

#include "Generator.hpp"
#include "Bus.hpp"
#include "Line.hpp"
#include "misc.hpp"

using namespace std;

class PowSys {

public:
    PowSys ();
    bool readData (string inputDir);

    int numGen;
    int numBus;
    int numLine;
    
    vector<Generator>   generators;
    vector<Bus>         buses;
    vector<Line>        lines;

private:
    
    bool readGeneratorData (string &inputDir);
    bool readBusData (string &inputDir);
    bool readLineData (string &inputDir);
    void postprocessing ();

    // helpers
    map<string, int> mapBusNameToIndex;    
};

#endif /* PowSys_hpp */
