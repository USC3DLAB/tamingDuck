//
//  Load.hpp
//  tamingDuck
//
//  Created by Semih Atakan on 12/17/17.
//  Copyright © 2017 Semih Atakan. All rights reserved.
//

#ifndef Load_hpp
#define Load_hpp

#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "misc.hpp"
#include "config.hpp"

using namespace std;

class Load {

public:
	Load ();
	bool read (string filepath);
	
	vector< vector<double> > regionalLoad;	// FirstIndex = regionId, SecondIndex = LoadForecast;	
};
#endif /* Load_hpp */
