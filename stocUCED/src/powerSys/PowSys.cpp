//
//  PowSys.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 12/14/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "PowSys.hpp"

extern runType runParam;

PowSys::PowSys () {}

bool PowSys::readData(string inputDir, string sysName) {

	bool status;

	name = sysName;
	path = inputDir + sysName + "/";

	// read generator data
	status = readGeneratorData(inputDir + sysName);
	if (!status) goto finalize;

	// read battery data
	status = readBatteryData(inputDir + sysName);
	
	// read bus data
	status = readBusData(inputDir + sysName);
	if (!status) goto finalize;

	// read line data
	status = readLineData(inputDir + sysName);
	if (!status) goto finalize;

	// postprocess
	postprocessing();

	finalize:
	if (status) {
		printf("Power system has been read successfully.\n");
	} else {
		printf("Error: Power system could not be read.\n");
	}
	return status;
}

bool PowSys::readBatteryData(string inputPath) {
	ifstream input;
	string temp_str;
	int batteryIndex = 0;
	
	/* Open file */
	bool status = open_file(input, inputPath + "/Batteries.csv");
	if (!status) {
		numBatteries = 0;
		printf("> Batteries.csv file not found (Optional).\n");
		return false;
	}
	
	// skip the headers
	safeGetline(input, temp_str);
	
	// read the data
	while (!input.eof()) {
		// create a generator
		Battery battery;
		battery.id = batteryIndex++;
		
		// skip the provided generator id
		move_cursor(input, delimiter);
		
		// generator name
		getline(input, battery.name, delimiter);
		
		// generator bus name
		getline(input, battery.connectedBusName, delimiter);
		
		// capacity
		input >> battery.maxCapacity;
		move_cursor(input, delimiter);
		
		// dissipation coefficient
		input >> battery.dissipationCoef;
		move_cursor(input, delimiter);
		
		// conversion losses
		input >> battery.chargingLossCoef;
		move_cursor(input, delimiter);
		input >> battery.dischargingLossCoef;		
		safeGetline(input, temp_str);
			
		// add the generator to the list
		batteries.push_back(battery);
	}
	input.close();
	
	numBatteries = (int)batteries.size();
	
	return true;
}

bool PowSys::readGeneratorData(string inputPath) {
	ifstream input;
	string temp_str;
	int genIndex = 0;

	/* Open file */
	bool status = open_file(input, inputPath + "/Generators.csv");
	if (!status)
		return false;

	// skip the headers
	safeGetline(input, temp_str);

	// read the data
	while (!input.eof()) {
		// create a generator
		Generator gen;
		gen.id = genIndex++;

		// skip the provided generator id
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

		// CO2 emissions
		input >> gen.CO2_emission_base;
		move_cursor(input, delimiter);

		input >> gen.CO2_emission_var;
		move_cursor(input, delimiter);

		// NOX emissions
		input >> gen.NOX_emission_base;
		move_cursor(input, delimiter);

		input >> gen.NOX_emission_var;
		move_cursor(input, delimiter);

		// SO2 emissions
		input >> gen.SO2_emission_base;
		move_cursor(input, delimiter);

		input >> gen.SO2_emission_var;
		move_cursor(input, delimiter);

		// must-run?
		input >> gen.isMustRun;
		move_cursor(input, delimiter);

		// must-use?
		input >> gen.isMustUse;
		move_cursor(input, delimiter);

		// day-ahead generator?
		input >> gen.isDAUCGen;
		safeGetline(input, temp_str);

		// add the generator to the list
		generators.push_back(gen);
	}
	input.close();

	numGen = (int)generators.size();

	return true;
}

bool PowSys::readBusData(string inputPath) {

	// open file
	ifstream input;
	bool status = open_file(input, inputPath + "/Buses.csv");
	if (!status)  return false;

	string temp_str;

	// skip the headers
	safeGetline(input, temp_str);

	int busIndex = 0;

	// read data
	while (!input.eof()) {
		// create a bus
		Bus bus;
		bus.id = busIndex;

		// bus name
		getline(input, bus.name, delimiter);

		// region id (may not be starting from 0)
		input >> bus.regionId;
		move_cursor(input, delimiter);

		// load percentage
		input >> bus.loadPercentage;
		safeGetline(input, temp_str);	// finalize reading

		// add bus to the list
		buses.push_back(bus);

		// add the index of the bus to the list
		mapBusNameToIndex.insert( pair<string, int> (bus.name, busIndex++) );
	}
	input.close();

	numBus = (int)buses.size();

	return true;
}

bool PowSys::readLineData(string inputPath) {

	// open file
	ifstream input;

	bool status = open_file(input, inputPath + "/Lines.csv");
	if (!status)  return false;

	string temp_str;

	// skip the headers
	safeGetline(input, temp_str);

	int lineIndex = 0;

	// read the data
	while (!input.eof()) {
		// create a line
		Line line;
		line.id = lineIndex++;

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
		safeGetline(input, temp_str);

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
 * Sets the connectedBus field of each battery
 * Sets the connectedBatteries field of each bus
 * Sets the connectedGenerators field of each bus
 * Scale solar/wind generator outputs
 * Scale generators' ramping capabilities
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
		
		// boost solar/wind generator outputs
		if (gen_ptr->type == Generator::SOLAR || gen_ptr->type == Generator::WIND) {
			gen_ptr->maxCapacity *= runParam.renewableCoef;
		}
		
		// adjust ramping coefficients
		if (gen_ptr->type != Generator::SOLAR && gen_ptr->type != Generator::WIND) {
			gen_ptr->rampUpLim *= runParam.rampingCoef;
			gen_ptr->rampDownLim *= runParam.rampingCoef;
		}
	}
	
	for (int bt=0; bt<numBatteries; bt++) {
		// get a handle of the generator
		Battery *bat_ptr = &(batteries[bt]);
		
		// get a handle of the bus
		Bus *bus_ptr = &(buses[ mapBusNameToIndex[bat_ptr->connectedBusName] ]);
		
		// set the connected bus of the generator
		bat_ptr->connectedBus = bus_ptr;
		
		// add the generator to the list of connected generators of the bus
		bus_ptr->connectedBatteries.push_back( bat_ptr );
	}
}
