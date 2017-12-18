//
//  PowSys.hpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#ifndef PowSys_hpp
#define PowSys_hpp

#include "misc.hpp"

#include "Generator.hpp"
#include "Bus.hpp"
#include "Line.hpp"

using namespace std;

class PowSys {

public:
    PowSys ();
    bool readData (string inputDir, string sysName);

    int numGen;
    int numBus;
    int numLine;
    
    vector<Generator>   generators;
    vector<Bus>         buses;
    vector<Line>        lines;

private:
    
    bool readGeneratorData (string inputPath);
    bool readBusData (string inputPath);
    bool readLineData (string inputPath);
    void postprocessing ();

    // helpers
    map<string, int> mapBusNameToIndex;    
};

#endif /* PowSys_hpp */
