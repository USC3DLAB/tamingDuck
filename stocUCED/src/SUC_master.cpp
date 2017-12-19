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
	/* TODO: Deal with the subproblem
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

    if (addCut)
	{
		// get the benders' cut
		double pi_b = me.sub.pi_b;
		vector<vector<double> >* pi_T = &me.sub.pi_T;
				
		IloExpr pi_Tx (me.env);
		for (int g=0; g<me.inst->G; g++) {
			for (int t=0; t<me.inst->T; t++) {
				if ( fabs((*pi_T)[g][t]) > 1e-10 ) {
					pi_Tx += (*pi_T)[g][t] * me.x[g][t];
				}
			}
		}
		
		IloRange BendersCut;
		if (isFeasible) {
			BendersCut = IloRange(me.env, pi_b, me.eta-pi_Tx);		// optimality cut
		} else {
			BendersCut = IloRange(me.env, pi_b, -pi_Tx);			// feasibility cut
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
	 */
}

SUCmaster::SUCmaster () {
	model = IloModel(env);
	cplex = IloCplex(env);
}

SUCmaster::~SUCmaster() {
	env.end();
}

void SUCmaster::preprocessing() {
	
	// initialize basic parameters
	numGen	   = inst->powSys->numGen;
	numLine    = inst->powSys->numLine;
	numBus     = inst->powSys->numBus;
	
	numPeriods	 = (probType == DayAhead) ? runParam.DA_numPeriods : runParam.ST_numPeriods;
	periodLength = (probType == DayAhead) ? runParam.DA_resolution : runParam.ST_resolution;
	
	numBaseTimePerPeriod = (int)round(periodLength) / runParam.ED_resolution;
}


/*
void master::formulate (instance &inst, int model_id) {
	
	// initialize the instance
	this->inst = &inst;
	
	// create the variables
	eta = IloNumVar (env);	// second-stage expectation approximation
    
    s = IloArray< IloNumVarArray > (env, inst.G);
	x = IloArray< IloNumVarArray > (env, inst.G);
	z = IloArray< IloNumVarArray > (env, inst.G);
	for (int g=0; g<inst.G; g++) {
		s[g] = IloNumVarArray(env, inst.T, 0, 1, ILOBOOL);
		x[g] = IloNumVarArray(env, inst.T, 0, 1, ILOBOOL);
		z[g] = IloNumVarArray(env, inst.T, 0, 1, ILOBOOL);
		
		char buffer[30];
		
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
    for (int g=0; g<inst.G; g++) {
        
        // t=0: generators are assumed to be turned on
        model.add( x[g][0] - 1 == s[g][0] - z[g][0] );
        
        // t>0
        for (int t=1; t<inst.T; t++) {
            model.add( x[g][t] - x[g][t-1] == s[g][t] - z[g][t] );
        }
    }
    
    // minimum uptime/downtime constraints
    for (int g=0; g<inst.G; g++)
    {
        // turn on inequalities
        for (int t=1; t<=inst.T; t++)
        {
            IloExpr lhs (env);
            for (int i = t-inst.min_u_time[g]+1; i<=t; i++)
            {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += 0;	// otherwise, generator is assumed to be operational (but turned on way earlier in the past)
            }
            model.add( lhs <= x[g][t-1] );
            lhs.end();
        }
    }
    
    for (int g=0; g<inst.G; g++)
    {
        // turn off inequalities
        for (int t=1; t<=inst.T; t++)
        {
            IloExpr lhs (env);
            for (int i = t-inst.min_d_time[g]+1; i<=t; i++)  {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += 0;	// otherwise, generator is assumed to be operational (but turned on way earlier in the past)
            }
            
            if (t-inst.min_d_time[g]-1 >= 0)	model.add( lhs <= 1 - x[g][t-inst.min_d_time[g]-1] );
            else								model.add( lhs <= 1 - 1 );	// assumed to be remaining on for a long, long, while
            lhs.end();
        }
    }
    
    // demand-based valid inequality
    // the capacities of operational generators must exceed the system demand at any point
    for (int t=0; t<inst.T; t++) {
        IloExpr expr (env);
        for (int g=0; g<inst.G; g++) {
            expr += inst.capacity[g] * x[g][t];
        }
        model.add( expr >= inst.aggDemand[t] );
        expr.end();
    }
    
    // must-run units must be committed
    for (int g=0; g<inst.G; g++) {
        if (inst.must_run[g]) {
            for (int t=0; t<inst.T; t++) {
                x[g][t].setLB(1);
            }
        }
    }
    
    // create the objective function
    IloExpr obj_func (env);
    
    for (int g=0; g<inst.G; g++) {
        for (int t=0; t<inst.T; t++) {
            obj_func += inst.start_cost[g] * s[g][t];						// start up cost
            obj_func += inst.no_load_cost[g] * x[g][t];                     // no-load cost
        }
    }
    obj_func += eta;
    
    // set the objective function
	model.add( IloMinimize(env, obj_func) );
	
    // formulate the subproblem
    sub.formulate(inst, model_id);

    /*************************************************************************
    /* DISABLED: Didn't improve the progress of the algorithm. Donno why..
     
    // compute a lower bound on eta
    double eta_lb = sub.computeLowerBound();
    eta.setLB( eta_lb );
    /************************************************************************
     
    // prepare the solver
	cplex.extract(model);

    /*************************************************************************
    /* DISABLED: Didn't improve the progress much.. reduced it in some cases.
    
    // assign priorities to state variables
    for (int g=0; g<inst.G; g++) {
        for (int t=0; t<inst.T; t++) {
            cplex.setPriority(x[g][t], 1);
        }
    }
    /************************************************************************
    
	cplex.use( LazySepCallback(env, *this) );
	cplex.setParam(IloCplex::Threads, 1);
    cplex.setParam(IloCplex::FPHeur, 2);
    cplex.setParam(IloCplex::HeurFreq, 10);
    cplex.setParam(IloCplex::LBHeur, 1);    // ??
    cplex.setParam(IloCplex::EpGap, 1e-2);
    cplex.setParam(IloCplex::Reduce, 0);    // due to callbacks
    cplex.setParam(IloCplex::PreLinear, 0); // due to callbacks
}
 */

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
void master::displayCommitments()
{
	try {
		for (int g=0; g<inst->G; g++) {
			for (int t=0; t<inst->T; t++) {
				cout << setprecision(0) << fixed << cplex.getValue(s[g][t] + x[g][t]) << " ";
			}
			cout << endl;
		}
	}
	catch (IloException &e) {
		cout << e << endl;
	}
}

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
