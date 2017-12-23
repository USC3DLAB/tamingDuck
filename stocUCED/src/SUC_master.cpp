//
//  SUC_master.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//

#include "SUC_master.hpp"

extern runType runParam;

IloCplex::Callback LazySepCallback(IloEnv env, SUCmaster & me) {
	return (IloCplex::Callback(new(env) SUCmaster::LazySepCallbackI(env, me)));
}

void SUCmaster::LazySepCallbackI::main()
{
	// get the solution
	vector< vector<bool> > state;
	resize_matrix(state, me.numGen, me.numPeriods);

    IloNumArray x_vals (me.env);
	for (int g=0; g<me.numGen; g++) {
		getValues(x_vals, me.x[g]);
		for (int t=0; t<me.numPeriods; t++) {
			(x_vals[t] >= 0.5) ? state[g][t] = true : state[g][t] = false;
		}
	}
    x_vals.end();

	// set the master solution in the subproblems
	me.sub.setMasterSoln(state);
	
	// solve the subproblems
	bool isFeasible = me.sub.solve();
	
    // check if we need to add a cut (feasible and not within optimality tolerances, or infeasible)
    bool addCut = true;
    if (isFeasible) {
        double abs_diff = fabs(me.sub.getRecourseObjValue() - getValue(me.eta));
        if ( abs_diff/(fabs(me.sub.getRecourseObjValue())+1e-14) < me.cplex.getParam(IloCplex::EpGap)
            || abs_diff < me.cplex.getParam(IloCplex::EpAGap) ) {
            addCut = false;
        }
    }

    if (addCut) {
		// get the Benders's cut coefs
		SUCsubprob::BendersCutCoefs *cutCoefs = &(me.sub.cutCoefs);
		
		// create the cut
		IloExpr pi_Tx (me.env);
		for (int g=0; g<me.numGen; g++) {
			for (int t=0; t<me.numPeriods; t++) {
				if ( fabs(cutCoefs->pi_T[g][t]) > 1e-10 ) {
					pi_Tx += cutCoefs->pi_T[g][t] * me.x[g][t];
				}
			}
		}
		
		IloRange BendersCut;
		if (isFeasible) {
			BendersCut = IloRange(me.env, cutCoefs->pi_b, me.eta-pi_Tx);	// optimality cut
		} else {
			BendersCut = IloRange(me.env, cutCoefs->pi_b, -pi_Tx);			// feasibility cut
		}
		
		// add the cut
		try {
			add(BendersCut).end();
		}
		catch (IloException &e) {
			cout << "Exception: " << e << endl;
			BendersCut.end();
		}
		pi_Tx.end();
	}
}

SUCmaster::SUCmaster () {
	model = IloModel(env);
	cplex = IloCplex(env);
}

SUCmaster::~SUCmaster() {
	env.end();
}

/****************************************************************************
 * preprocessing
 * - Initializes basic parameters
 * - Fills the following containers, according to model units & assumptions:
 *   minUpTimePeriods
 *   minDownTimePeriods
 *   expCapacity
 *   aggLoads
 * - If it is a DA-UC problem, sets the capacity and minGenerationReq to 0
 * for generators which must be scheduled at ST-UC.
 *
 * - Assumption 1: Min generation requirement must be less than ramping
 * rates, otherwise, generators cannot switch on/off.
 * - Assumption 2: Min up/down times must at least be 1 period.
 ****************************************************************************/
void SUCmaster::preprocessing ()
{
	/* basic parameters */
	numGen	   = inst->powSys->numGen;
	numLine    = inst->powSys->numLine;
	numBus     = inst->powSys->numBus;
	
	numPeriods	 = (probType == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
	periodLength = (probType == DayAhead) ? runParam.DA_resolution : runParam.ST_resolution;
	
	numBaseTimePerPeriod = (int)round(periodLength) / runParam.ED_resolution;
	
	
	/* initialize containers */
	minUpTimePeriods.resize(numGen);					// minimum uptime in periods
	minDownTimePeriods.resize(numGen);					// minimum downtime in periods
	resize_matrix(expCapacity, numGen, numPeriods);		// mean generator capacities
	sysLoad.resize(numPeriods);							// aggregated system load
	
	/* Min Up/Down */
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		minUpTimePeriods[g]		= round(genPtr->minUpTime * 60.0/periodLength);
		minDownTimePeriods[g]	= round(genPtr->minDownTime * 60.0/periodLength);
		
		if (minUpTimePeriods[g] < 1)	minUpTimePeriods[g] = 1;
		if (minDownTimePeriods[g] < 1)	minDownTimePeriods[g] = 1;
	}
	
	/* Mean Generator Capacity */
	auto dataPtr = (probType == DayAhead) ? &(inst->stocObserv[0]) : &(inst->stocObserv[1]);
	
	for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst->powSys->generators[g]);
		
		auto it = dataPtr->mapVarNamesToIndex.find(genPtr->name);
		if ( it != dataPtr->mapVarNamesToIndex.end() ) {
			/* random supply */
			for (int t=0; t<numPeriods; t++) {
				expCapacity[g][t] = dataPtr->vals[0][t*numBaseTimePerPeriod][it->second];	// in MWs
			}
		} else {
			/* deterministic supply */
			fill(expCapacity[g].begin(), expCapacity[g].end(), genPtr->maxCapacity);	// in MWs
		}
	}
	
	/* Misc: Generators which are scheduled in ST-UC, must have 0 min-generation-requirement & capacity in DA-UC problem */
	if (probType == DayAhead) {
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst->powSys->generators[g]);
			
			if (!genPtr->isBaseLoadGen) {
				fill(expCapacity[g].begin(), expCapacity[g].end(), 0.0);
			}
		}
	}
	
	/* Load */
	dataPtr = (probType == DayAhead) ? &(inst->detObserv[0]) : &(inst->detObserv[1]);
	fill( sysLoad.begin(), sysLoad.end(), 0.0 );		// initialize to 0
	for (int t=0; t<numPeriods; t++) {
		for (int r=0; r<dataPtr->numVars; r++) {
			sysLoad[t] += dataPtr->vals[0][t*numBaseTimePerPeriod][r];
		}
	}
}


void SUCmaster::formulate (instance &inst, ProblemType probType, ModelType modelType, int beginMin) {
	
	/* Initializations */
	this->inst		= &inst;
	this->beginMin	= beginMin;
	this->probType	= probType;
	this->modelType = modelType;
	
	/* Prepare Model-Dependent Input Data */
	preprocessing();

	/* Master Formulation */
	// create the variables
	eta = IloNumVar (env);	// second-stage expectation approximation
    
    s = IloArray< IloNumVarArray > (env, numGen);
	x = IloArray< IloNumVarArray > (env, numGen);
	z = IloArray< IloNumVarArray > (env, numGen);
	for (int g=0; g<numGen; g++) {
		s[g] = IloNumVarArray(env, numPeriods, 0, 1, ILOBOOL);
		x[g] = IloNumVarArray(env, numPeriods, 0, 1, ILOBOOL);
		z[g] = IloNumVarArray(env, numPeriods, 0, 1, ILOBOOL);
		
		sprintf(buffer, "s_%d", g);
		s[g].setNames(buffer);
		sprintf(buffer, "x_%d", g);
		x[g].setNames(buffer);
		sprintf(buffer, "z_%d", g);
		z[g].setNames(buffer);
		
		model.add(s[g]);
		model.add(x[g]);
		model.add(z[g]);
	}
	
	// create the constraints

    // state constraints
    for (int g=0; g<numGen; g++) {
        
        // t=0: generators are assumed to be turned on
        model.add( x[g][0] - getGenState(g,-1) == s[g][0] - z[g][0] );
        
        // t>0
        for (int t=1; t<numPeriods; t++) {
            model.add( x[g][t] - x[g][t-1] == s[g][t] - z[g][t] );
        }
    }
    
    // minimum uptime/downtime constraints
    for (int g=0; g<numGen; g++)
    {
        // turn on inequalities
        for (int t=1; t<=numPeriods; t++)
        {
            IloExpr lhs (env);
            for (int i = t-minUpTimePeriods[g]+1; i<=t; i++)
            {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += max(0, getGenState(g,i-1) - getGenState(g,i-2));
            }
            model.add( lhs <= x[g][t-1] );
            lhs.end();
        }
    }
    
    for (int g=0; g<numGen; g++)
    {
        // turn off inequalities
        for (int t=1; t<=numPeriods; t++)
        {
            IloExpr lhs (env);
            for (int i = t-minDownTimePeriods[g]+1; i<=t; i++)  {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += max(0, getGenState(g,i-1) - getGenState(g,i-2));
            }
            
			if (t-minDownTimePeriods[g]-1 >= 0)	model.add( lhs <= 1 - x[g][t-minDownTimePeriods[g]-1] );
			else								model.add( lhs <= 1 - getGenState(g, (t-minDownTimePeriods[g]-1)));

			lhs.end();
        }
    }
    
    // demand-based valid inequality
	// TODO: Revisit this 'valid' inequality
    // the capacities of operational generators must exceed the system demand at any point
/*    for (int t=0; t<numPeriods; t++) {
        IloExpr expr (env);
        for (int g=0; g<numGen; g++) {
            expr += expCapacity[g] * x[g][t];
        }
        model.add( expr >= inst.aggDemand[t] );
        expr.end();
    }
 */
    
	// must-run units must be committed
	for (int g=0; g<numGen; g++) {
		if ( inst.powSys->generators[g].isMustRun ) {
			for (int t=0; t<numPeriods; t++) {
				x[g][t].setBounds(1, 1);
			}
		}
	}
	
	// commit the right set of generators for the right type of problem
	if (probType == DayAhead) {
		/* ST-UC generators will not produce in the DA-UC problem. */
		for (int g=0; g<numGen; g++) {
			if (!(inst.powSys->generators[g].isBaseLoadGen)) {
				for (int t=0; t<numPeriods; t++) {
					x[g][t].setBounds(0, 0);
				}
			}
		}
	}
	else {
		/* All generators will produce in the ST-UC problem. Commitment
		 * decisions of DA-UC generators will be read from the solution. */
		double genState;
		for (int g=0; g<numGen; g++) {
			Generator *genPtr = &(inst.powSys->generators[g]);
			
			if (genPtr->isBaseLoadGen) {
				for (int t=0; t<numPeriods; t++) {
					genState = getGenState(g, t);
					x[g][t].setBounds(genState, genState);
				}
			}
		}
	}
	
    // create the objective function
    IloExpr obj (env);
    
    for (int g=0; g<numGen; g++) {
		Generator *genPtr = &(inst.powSys->generators[g]);
		
        for (int t=0; t<numPeriods; t++) {
            obj += genPtr->startupCost * s[g][t];					// start up cost
            obj += genPtr->noLoadCost*periodLength/60.0 * x[g][t];	// no-load cost
        }
    }
    obj += eta;
    
    // set the objective function
	model.add( IloMinimize(env, obj) );
	
    // formulate the subproblem
	sub.formulate(inst, probType, modelType, beginMin);

    /*************************************************************************
     * DISABLED: Didn't improve the progress of the algorithm. Donno why..
     *
     * // compute a lower bound on eta
     * double eta_lb = sub.computeLowerBound();
     * eta.setLB( eta_lb );
     ************************************************************************/
     
    // prepare the solver
	cplex.extract(model);

    /*************************************************************************
     * DISABLED: Didn't improve the progress much.. reduced it in some cases.
     *
     * // assign priorities to state variables
     * for (int g=0; g<inst.G; g++) {
	 *    for (int t=0; t<inst.T; t++) {
     *         cplex.setPriority(x[g][t], 1);
     *     }
     * }
     ************************************************************************/
    
	cplex.use( LazySepCallback(env, *this) );
	cplex.setParam(IloCplex::Threads, 1);
    cplex.setParam(IloCplex::FPHeur, 2);
    cplex.setParam(IloCplex::HeurFreq, 10);
    cplex.setParam(IloCplex::LBHeur, 1);    // ??
    cplex.setParam(IloCplex::EpGap, 1e-2);
    cplex.setParam(IloCplex::Reduce, 0);    // due to callbacks
    cplex.setParam(IloCplex::PreLinear, 0); // due to callbacks
}

bool SUCmaster::solve () {
	try{
        bool status = cplex.solve();
        
		cout << "Optimization is completed with status " << cplex.getCplexStatus() << endl;
        cout << "Obj = \t" << cplex.getObjValue() << endl;
        cout << "LB = \t" << cplex.getBestObjValue() << endl;
        
		return status;
	}
	catch (IloException &e) {
		cout << e << endl;
		return false;
	}
}

/*
void master::warmStart(Solution &soln)
{
    this->soln = soln;
    
    IloNumVarArray vars (env);
    IloNumArray vals (env);
    
    for (int g=0; g<inst->G; g++) {
        for (int t=0; t<inst->T; t++) {
            vars.add( s[g][t] );
            vals.add( soln.s[g][t] );
            
            vars.add( x[g][t] );
            vals.add( soln.x[g][t] );

            vars.add( z[g][t] );
            vals.add( soln.z[g][t] );
        }
    }
    
    cplex.addMIPStart(vars, vals);
    
    vars.end();
    vals.end();
}

void master::warmStart(vector<Solution> &solnPool)
{
    this->soln = solnPool[0];
    
    IloNumVarArray vars (env);
    IloNumArray vals (env);
    
    for (int sol=0; sol<solnPool.size(); sol++)
    {
        for (int g=0; g<inst->G; g++) {
            for (int t=0; t<inst->T; t++) {
                vars.add( s[g][t] );
                vals.add( solnPool[sol].s[g][t] );
                
                vars.add( x[g][t] );
                vals.add( solnPool[sol].x[g][t] );
                
                vars.add( z[g][t] );
                vals.add( solnPool[sol].z[g][t] );
            }
        }
        cplex.addMIPStart(vars, vals);
    }
    
    vars.end();
    vals.end();
}

*/

/****************************************************************************
 * getGenState
 * - Converts the model period, into the desired component of the Solution
 * object. Returns the recorded status of the generator.
 ****************************************************************************/
bool SUCmaster::getGenState(int genId, int period) {
	// which Solution component is requested?
	int reqSolnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// return the requested generator state
	if (reqSolnComp < 0) {										// all generators are assumed to be online, for a long time, at t=0.
		return true;
	}
	else if (reqSolnComp < inst->solution.x[genId].size()) {	// return the corresponding solution
		return round(inst->solution.x[genId][reqSolnComp]);
	}
	else {														// error
		cout << "Error: You cannot access a solution component that is beyond the planning horizon" << endl;
		exit(-1);
	}
}

/****************************************************************************
 * setGenState
 * - Fills the (genId, correspondingComponent) of the Solution.x object.
 ****************************************************************************/
void SUCmaster::setGenState(int genId, int period, double value) {
	// which Solution component is being set?
	int solnComp = beginMin/runParam.ED_resolution + period*numBaseTimePerPeriod;
	
	// set the solution
	if (solnComp >= 0 && solnComp < inst->solution.x[genId].size()) {
		for (int t=solnComp; t<solnComp+numBaseTimePerPeriod; t++) {
			inst->solution.x[genId][t] = value;
		}
	}
	else {
		cout << "Error: Setting generator state out of the bounds of the horizon" << endl;
		exit(-1);
	}
}
