//
//  PowSys.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "PowSys.hpp"

PowSys::PowSys () {}

bool PowSys::readData(string inputDir) {
    
    bool status;
    
    // read generator data
    status = readGeneratorData(inputDir);
    if (!status) return false;
    
    // read bus data
    status = readBusData(inputDir);
    if (!status) return false;
    
    // read line data
    status = readLineData(inputDir);
    if (!status) return false;
    
    // finalize
    postprocessing();
    
    return true;
}

bool PowSys::readGeneratorData(string &inputDir) {
    // open file
    ifstream input;
    bool status = open_file(input, inputDir + "Generators.csv");
    if (!status)  return false;
    
    // read file
    
    // skip the headers
    move_cursor(input, eoline);
    
    while (!input.eof()) {
        // create a generator
        Generator gen;
        
        // get generator id
        input >> gen.id;
        move_cursor(input, delimiter);
        
        // generator name
        getline(input, gen.name, delimiter);
        
        // generator type
        string generatorTypeName;
        getline(input, generatorTypeName, delimiter);
        gen.setType(generatorTypeName);
        
        // generator bus name
        getline(input, gen.connectedBusName, delimiter);
        
        // capacity
        input >> gen.maxCapacity;
        move_cursor(input, delimiter);
        
        // variable cost
        input >> gen.variableCost;
        move_cursor(input, delimiter);
        
        // no-load cost
        input >> gen.noLoadCost;
        move_cursor(input, delimiter);
        
        // ramp down
        input >> gen.rampDownLim;
        move_cursor(input, delimiter);
        
        // ramp up
        input >> gen.rampUpLim;
        move_cursor(input, delimiter);
        
        // min down time
        input >> gen.minDownTime;
        move_cursor(input, delimiter);
        
        // min up time
        input >> gen.minUpTime;
        move_cursor(input, delimiter);
        
        // min stable time
        input >> gen.minGenerationReq;
        move_cursor(input, delimiter);
        
        // start up cost
        input >> gen.startupCost;
        move_cursor(input, delimiter);
        
        // must-run?
        input >> gen.isMustRun;
        move_cursor(input, delimiter);
        
        // day-ahead generator?
        input >> gen.isBaseLoadGen;
        move_cursor(input, eoline);
        
        // add the generator to the list
        generators.push_back(gen);
    }
    input.close();
    
    numGen = (int)generators.size();
    
    return true;
}

bool PowSys::readBusData(string &inputDir) {
  
    // open file
    ifstream input;
    bool status = open_file(input, inputDir + "Buses.csv");
    if (!status)  return false;
    
    // skip the headers
    move_cursor(input, '\n');
    //TODO: These endline tokens will be future cause of problem
    
    int busIndex = 0;
    
    // read data
    while (!input.eof()) {
        // create a bus
        Bus bus;
        
        // bus name
        getline(input, bus.name, delimiter);

        // region id (may not be starting from 0)
        input >> bus.regionId;
        move_cursor(input, delimiter);
       
        // load percentage
        input >> bus.loadPercentage;
        move_cursor(input, '\n');
        
        // add bus to the list
        buses.push_back(bus);
        
        // add the index of the bus to the list
        mapBusNameToIndex.insert( pair<string, int> (bus.name, busIndex++) );
    }
    input.close();
    
    numBus = (int)buses.size();

    return true;
}

bool PowSys::readLineData(string &inputDir) {
    
    // open file
    ifstream input;

    bool status = open_file(input, inputDir + "Lines.csv");
    if (!status)  return false;
    
    // skip the headers
    move_cursor(input, eoline);
    
    // read the data
    while (!input.eof()) {
        // create a line
        Line line;
        
        // name
        getline(input, line.name, delimiter);
        
        // origin and destination buses
        string busName;
        
        getline(input, busName, delimiter);
        line.orig = &(buses[mapBusNameToIndex[busName]]);
        
        getline(input, busName, delimiter);
        line.dest = &(buses[mapBusNameToIndex[busName]]);
        
        // max flow
        input >> line.maxFlowLim;
        move_cursor(input, delimiter);
        
        // min flow
        input >> line.minFlowLim;
        move_cursor(input, delimiter);
                
        // read susceptance
        input >> line.susceptance;
        move_cursor(input, eoline);
        
        // add the line to the list
        lines.push_back(line);
    }
    input.close();
    
    numLine = (int) lines.size();

    return true;
}

/****************************************************************************
 * postprocessing
 * Sets the connectedBus field of each generator
 * Sets the connectedGenerators field of each bus
 ****************************************************************************/
void PowSys::postprocessing() {
    
    for (int g=0; g<numGen; g++) {
        // get a handle of the generator
        Generator *gen_ptr = &(generators[g]);
        
        // get a handle of the bus
        Bus *bus_ptr = &(buses[ mapBusNameToIndex[gen_ptr->connectedBusName] ]);
        
        // set the connected bus of the generator
        gen_ptr->connectedBus = bus_ptr;
        
        // add the generator to the list of connected generators of the bus
        bus_ptr->connectedGenerators.push_back( gen_ptr );
    }
}
