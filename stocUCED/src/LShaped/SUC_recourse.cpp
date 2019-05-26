///
//  SUC_recourse.cpp
//  tamingDuck
//
//  Created by Semih Atakan on 2/17/18.
//  Copyright © 2018 Semih Atakan. All rights reserved.
//

#include "SUC_recourse.hpp"

extern runType runParam;

SUCrecourse::SUCrecourse () {
	subprobs.resize(numThreads);
}

SUCrecourse::~SUCrecourse() {}

void SUCrecourse::formulate(instance &inst, ProblemType probType, ModelType modelType, int beginMin, int rep, IloArray<IloNumArray> &masterSoln, vector<int> &rndPermutation, vector<vector<double>> &expCapacity) {
	numScen = runParam.numLSScen;
	
	int numGen = inst.powSys->numGen;
	int numPeriods = (probType == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
	int numBatteries = inst.powSys->numBatteries;

	/* formulate the subproblems */
	for (int k=0; k<subprobs.size(); k++) {
		subprobs[k].formulate(inst, probType, modelType, beginMin, rep, masterSoln, expCapacity);
	}
	
	// allocate mem
	objValues.resize(numScen);
	cutCoefs.resize(numScen);
	for (int s=0; s<numScen; s++) {
		cutCoefs[s].initialize(numGen, numPeriods);
	}
	
	resize_matrix(initGens, numScen, numGen);
	btStates.resize(numScen);
	for (int s=0; s<numScen; s++) {
		resize_matrix(btStates[s], numBatteries, numPeriods);
	}
	
	this->rndPermutation = rndPermutation;
}

/****************************************************************************
 * solve
 * - Solves the subproblems, gets the cut coefficients, obj values, and 
 * initial generation amounts.
 * - The function is implemented in parallel, if BOOST_PARALLEL_LIBS is defined
 ****************************************************************************/
#ifdef BOOST_PARALLEL_LIBS
bool SUCrecourse::solve() {
	bool status = true;
	
	// shuffle subproblem solving sequence
	vector<int> randomized_order (numScen);
	for (int s=0; s<randomized_order.size(); s++) randomized_order[s] = s;
	random_shuffle(randomized_order.begin(), randomized_order.end());
	
	/***** Parallel programming stuff (START) *****/
	boost::asio::io_service io_service;				// create an io_service
	boost::asio::io_service::work work(io_service);	// and some work to stop its run() function from exiting if it has nothing else to do
	boost::thread_group threads;					// start some worker threads
	
	thread_map.clear();								// create new threads, get their ids into a map
	for (int k=0; k<numThreads; k++)	{
		boost::thread *new_thread = threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
		thread_map.insert( pair<boost::thread::id, unsigned short> (new_thread->get_id(), k) );
	}
	/***** Parallel programming stuff (END)   *****/
	
	// solve subproblems
	for (int s=0; s<numScen; s++) {
		io_service.post( boost::bind(&SUCrecourse::solveOneSubproblem, boost::ref(*this),
									 rndPermutation[ randomized_order[s] ],
									 boost::ref(cutCoefs[ randomized_order[s] ]),
									 boost::ref(objValues[ randomized_order[s] ]),
									 boost::ref(initGens[ randomized_order[s] ]),
									 boost::ref(btStates[ randomized_order[s] ])) );
	}

	/***** Parallel programming stuff (START) *****/
	work.~work();		// let the io_service shutdown by removing the thread
	threads.join_all();
	/***** Parallel programming stuff (END)   *****/
	
	/* check feasibility */
	for (int s=0; s<numScen; s++) {
		if (objValues[s] == INFINITY)	{
			status = false;
			break;
		}
	}

	return status;
}

void SUCrecourse::solveOneSubproblem (int s, BendersCutCoefs &cutCoefs, double &objValue, vector<double> &initGens, vector<vector<double>> &btStates) {
	subprobs[ thread_map[boost::this_thread::get_id()] ].solve(s, cutCoefs, objValue, initGens, btStates);
}

#else
bool SUCrecourse::solve() {
	bool status = true;
	
	// process second-stage scenario subproblems
	for (int s=0; s<numScen; s++) {
		status = subprobs[0].solve(rndPermutation[s], cutCoefs[s], objValues[s], initGens[s], btStates[s]);
		
		if (!status) {
			break;
		}
	}
	
	return status;
}
#endif

bool SUCrecourse::solveMeanProb() {
	return subprobs[0].solve(-1, cutCoefs[0], objValues[0], initGens[0], btStates[0]);
}

/****************************************************************************
 * setMasterSoln
 * - Sets the master solution at each subprob math model (or solver) object
 ****************************************************************************/
void SUCrecourse::setMasterSoln() {
	for (int k=0; k<subprobs.size(); k++) {
		subprobs[k].setMasterSoln();
	}
}

/****************************************************************************
 * getObjValue
 * - Returns the recourse problem's objective value (avg. of subproblem objs)
 ****************************************************************************/
double SUCrecourse::getObjValue() {
	double objVal = 0;
	for (int s=0; s<numScen; s++) {
		objVal += (1.0/double(numScen)) * objValues[s];
	}
	return objVal;
}

/****************************************************************************
 * getScenObjValue
 * - Returns the objective value of s^th subproblem
 ****************************************************************************/
double SUCrecourse::getScenObjValue(int s) {
	return objValues[s];
}

/****************************************************************************
 * getExpInitGen
 * - Computes the sample-average of period-0 generation levels, which were
 * recorded in the Benders' subproblems.
 ****************************************************************************/
vector<double> SUCrecourse::getExpInitGen() {
	vector<double> expInitGen (initGens[0].size(), 0.0);	// # of generators-many
	
	for (int s=0; s<numScen; s++) {
		for (int g=0; g<initGens[0].size(); g++) {
			expInitGen[g] += 1.0/(double)numScen * initGens[s][g];
		}
	}
	return expInitGen;
}

/****************************************************************************
 * getExpBtState
 * - Computes the sample-average of battery state levels, which were
 * recorded in the Benders' subproblems.
 ****************************************************************************/
vector<vector<double>> SUCrecourse::getExpBtState() {
	vector<vector<double>> expBtState (btStates[0].size(), vector<double> (btStates[0][0].size(), 0.0));	// # of (batteries x periods)-many
	
	for (int s=0; s<numScen; s++) {
		for (int bt=0; bt<btStates[0].size(); bt++) {
			for (int t=0; t<btStates[0][0].size(); t++) {
				expBtState[bt][t] += 1.0/(double)numScen * btStates[s][bt][t];
			}
		}
	}
	return expBtState;
}

/****************************************************************************
 * getInfeasScenIndex
 * - Returns the index of the first infeasible-scenario encountered.
 ****************************************************************************/
int SUCrecourse::getInfeasScenIndex() {
	for (int s=0; s<numScen; s++) {
		if (objValues[s] == INFINITY) {
			return s;
		}
	}
	return -1;
}
