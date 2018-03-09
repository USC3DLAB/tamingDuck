//
//  SUC_recourse.hpp
//  tamingDuck
//
//  Created by Semih Atakan on 2/17/18.
//  Copyright © 2018 Semih Atakan. All rights reserved.
//

#ifndef SUC_recourse_hpp
#define SUC_recourse_hpp

#include <iostream>
#include "SUC_subprob.hpp"
#include "config.hpp"		// must have this before defining boost libraries

#ifdef BOOST_PARALLEL_LIBS	// the boost library is added for multi-threading
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/asio.hpp>
#endif

class SUCrecourse {
	friend class SUCmaster;
	friend class LShapedCallback;
public:
	SUCrecourse ();
	~SUCrecourse();
	
	void formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep, IloArray<IloNumArray> &masterSoln, vector<int> &rndPermutation, vector<vector<double>> &expCapacity);
	bool solve ();
	bool solveMeanProb ();
	
	void setMasterSoln ();
	
	int	getInfeasScenIndex ();
	
	double getObjValue();
	double getScenObjValue(int s);
	
	vector<double> getExpInitGen ();
private:
	/* Benders' cuts */
	vector<BendersCutCoefs> cutCoefs;

	/* Subproblems */
	vector<SUCsubprob> subprobs;
	
	/* Scenario-wise Containers */
	vector<double> objValues;
	vector< vector<double> > initGens;
	
	/* Parallelization */
	const int numThreads = LShapedSubprobThreads;
#ifdef BOOST_PARALLEL_LIBS
	map<boost::thread::id, unsigned short> thread_map;
	void solveOneSubproblem (int s, BendersCutCoefs &cutCoefs, double &objValue, vector<double> &initGens);
#endif
	
	/* Misc */
	int numScen;	
	int infeasScenIdx;
	
	vector<int> rndPermutation;
};

#endif /* SUC_recourse_hpp */
