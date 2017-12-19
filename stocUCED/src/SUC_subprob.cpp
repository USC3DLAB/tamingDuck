//
//  subprob.cpp
//  Stoch-UC
//
//  Created by Semih Atakan on 10/31/17.
//  Copyright © 2017 University of Southern California. All rights reserved.
//
/*
#include "subprob.hpp"

subprob::subprob () {
	timer_t = 0.0;
	
	model = IloModel (env);
	cplex = IloCplex (env);
	
	/* TODO: Commented by HG after encountering "Invalid arguments candidates are: void setOut(? &)
	cplex.setOut(env.getNullStream()); *
	cplex.setParam(IloCplex::PreInd, 0);
	cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
}

subprob::~subprob () {
	env.end();
}

void subprob::formulate (instance &inst, int model_id)
{
	// get a handle on the instance and model type
	this->inst = &inst;
	this->model_id = model_id;
	
	// some initializations
	cons  = IloRangeArray (env);
	duals = IloNumArray (env);
	resize_matrix(pi_T, inst.G, inst.T);
		
	// formulate the desired problem
	switch (model_id) {
		case 0:
			formulate_aggregate_system();
			break;
		case 1:
			formulate_nodebased_system();
			break;
		default:
			std::cout << "Invalid model ID" << std::endl;
			break;
	}
}

void subprob::formulate_production()
{
	/** Decision Variables **
	p	  = IloArray< IloNumVarArray > (env, inst->G);	// production amounts
	for (int g=0; g<inst->G; g++) {
		p[g]	 = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
		
		char buffer[30];
		sprintf(buffer, "p_%d", g);
		p[g].setNames(buffer);
	}

	/** Constraints **
 
	// capacity constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++) {
			IloRange con;
			if (inst->must_run[g]) {
                con = IloRange (env, inst->capacity[g], p[g][t], inst->capacity[g]);
			} else {
                con = IloRange (env, -IloInfinity, p[g][t], inst->capacity[g]);
			}
			cons.add( con );
		}
	}
	
	// minimum production constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++) {
			IloRange con (env, inst->min_prod_lim[g], p[g][t]);
			cons.add( con );
		}
	}
	
	// ramp up constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=1; t<inst->T; t++) {
			IloRange con (env, -IloInfinity, p[g][t] - p[g][t-1], inst->ramp_u_lim[g]);
			cons.add( con );
		}
	}
	
	// ramp down constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=1; t<inst->T; t++) {
			IloRange con (env, -IloInfinity, p[g][t-1] - p[g][t], inst->ramp_d_lim[g]);
			cons.add( con );
		}
	}
}

void subprob::formulate_nodebased_system()
{
	formulate_production();
	
	/** Decision Variables **
    L = IloArray< IloNumVarArray > (env, inst->B);
	IloArray< IloNumVarArray > T (env, inst->B);		// phase angles
	for (int b=0; b<inst->B; b++) {
		L[b] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);  // load shedding
		T[b] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
		
		char buffer[30];
		
		sprintf(buffer, "L_%d", b);
		L[b].setNames(buffer);
		
		sprintf(buffer, "T_%d", b);
		T[b].setNames(buffer);
	}
	
	IloArray< IloNumVarArray > F (env, inst->L);		// flows
	for (int l=0; l<inst->L; l++) {
		F[l] = IloNumVarArray(env, inst->T, -IloInfinity, IloInfinity, ILOFLOAT);
		
		char buffer[30];
		
		sprintf(buffer, "F_%d", l);
		F[l].setNames(buffer);
	}
	
	/** Constraints **
	// Flow upper and lower bounds
	for (int l=0; l<inst->L; l++) {			// upper bound
		for (int t=0; t<inst->T; t++) {
			IloRange con(env, -IloInfinity, F[l][t], inst->max_arc_flow[l]);
			cons.add(con);
		}
	}
	for (int l=0; l<inst->L; l++) {			// lower bound
		for (int t=0; t<inst->T; t++) {
			IloRange con(env, inst->min_arc_flow[l], F[l][t]);
			cons.add(con);
		}
	}
	
	// Phase-angle upper bounds
	for (int b=0; b<inst->B; b++) {
		for (int t=0; t<inst->T; t++) {
			IloRange con(env, -IloInfinity, T[b][t], 2*inst->pi);
			cons.add(con);
		}
	}
	
	// DC-approximation to AC flow
	for (int l=0; l<inst->L; l++) {
		for (int t=0; t<inst->T; t++) {
			int orig = inst->arcs[l].first;
			int dest = inst->arcs[l].second;
			
			model.add( F[l][t] == inst->susceptance[l] * (T[orig][t] - T[dest][t]) );
		}
	}
	
	// Flow balance
	for (int b=0; b<inst->B; b++) {
		for (int t=0; t<inst->T; t++) {
			
			IloExpr expr (env);
			
			// production
			for (int g=0; g<inst->G; g++) {
				if (inst->generator_loc[g] == b)    expr += p[g][t];
			}

			// incoming arcs
			for (int l=0; l<inst->L; l++) {
				if (inst->arcs[l].second == b)      expr += F[l][t];
			}
			
			// outgoing arcs
			for (int l=0; l<inst->L; l++) {
				if (inst->arcs[l].first == b)       expr -= F[l][t];
			}
			
			// load shedding
			expr += L[b][t];
			
			// constraint
			IloRange con(env, inst->demand[b][t], expr, inst->demand[b][t]);
			cons.add(con);
			
			expr.end();
		}
	}
	
	
	/** Objective Function **
	IloExpr obj_func (env);
	
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++) {
			obj_func += inst->var_cost[g] * p[g][t];
		}
	}
	
	for (int b=0; b<inst->B; b++) {
		for (int t=0; t<inst->T; t++) {
			obj_func += inst->load_shedding_penalty * L[b][t];
		}
	}
	
	/** Finalize **
	model.add( IloMinimize(env, obj_func) );
	model.add( cons );
	cplex.extract(model);
	
}

void subprob::formulate_aggregate_system()
{
	formulate_production();

	/** Decision Variables **
	IloNumVarArray L (env, inst->T, 0, IloInfinity, ILOFLOAT);	// load shedding

	/** Constraints **
	// aggregated demand constraints
	for (int t=0; t<inst->T; t++) {
		IloExpr expr (env);
		for (int g=0; g<inst->G; g++)	expr += p[g][t];
		expr += L[t];
		IloRange con (env, inst->aggDemand[t], expr);
		cons.add( con );
	}
	
	/** Objective Function **
	IloExpr obj_func (env);
	
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++) {
			obj_func += inst->var_cost[g] * p[g][t];
		}
	}
	
	for (int t=0; t<inst->T; t++) {
		obj_func += inst->load_shedding_penalty * L[t];
	}
	
	// prepare the model and the solver
	model.add( IloMinimize(env, obj_func) );
	model.add( cons );
	cplex.extract(model);
}

void subprob::setMasterSoln ( vector< vector<bool> > & gen_stat )
{
	this->gen_stat = &gen_stat;
	
	int c=0;
	
	// capacity constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++, c++) {
            if (inst->must_run[g]) {
                cons[c].setBounds( gen_stat[g][t] * inst->capacity[g], gen_stat[g][t] * inst->capacity[g]);
            } else {
                cons[c].setUB( gen_stat[g][t] * inst->capacity[g] );
            }
		}
	}

	// production amounts
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++, c++) {
			cons[c].setLB( gen_stat[g][t] * inst->min_prod_lim[g] );
		}
	}
    
    // rest of the constraints are not a function of x
}

bool subprob::solve() {
	
	reset_cut_coefs();
	recourse_obj_val = 0;
	
	bool status = true;
	
	// process second-stage scenario subproblems
	for (int s=0; s<inst->S; s++)
	{
		setup_subproblem (s);			// prepare the subproblem
		
		status = cplex.solve();			// solve the subproblem
		if (!status) {
			
			if ( cplex.getCplexStatus() == IloCplex::Infeasible ) {
				get_feasibility_cut_coefs(s);
				recourse_obj_val = INFINITY;
			}
			else {
				cout << "!! Error: Cplex returned the status " << cplex.getCplexStatus() << endl;
			}
			
			return false;
		}
		
		recourse_obj_val += inst->sceProb[s] * cplex.getObjValue();	// get the obj value
		
		update_optimality_cut_coefs(s);	// update benders cut coefs
	}
	
	return true;
}

void subprob::setup_subproblem(int &s)
{
	// iterate over capacity constraints
	int c=0;
	
	// capacity constraints
	for (int g=0; g<inst->G; g++) {
		
		map<int, vector<vector<double> > >::iterator it = inst->rndSupply.find(g);
		
		if ( it != inst->rndSupply.end() ) {
			for (int t=0; t<inst->T; t++, c++) {	// Note: c is iterated in the secondary-loops
				if ((*gen_stat)[g][t]) {
                    if ( inst->must_run[g] ) {
                        cons[c].setBounds( inst->rndSupply[g][s][t], inst->rndSupply[g][s][t]);
                    } else {
                        cons[c].setUB ( inst->rndSupply[g][s][t] );
                    }
				}
			}
		}
		else {
			c += inst->T;
		}
	}    
}

void subprob::update_optimality_cut_coefs(int &s)
{
	// get dual multipliers
	cplex.getDuals(duals, cons);
	
	// iterate over the constraints to write the Benders' cut
	int c=0;
	
	// capacity constraints
	for (int g=0; g<inst->G; g++) {
		map<int, vector<vector<double> > >::iterator it = inst->rndSupply.find(g);
		
		if ( it != inst->rndSupply.end() ) {
			for (int t=0; t<inst->T; t++, c++) {	// Note: c is iterated in the secondary-loops
				pi_T[g][t] += inst->sceProb[s] * duals[c] * ( inst->rndSupply[g][s][t] );
			}
		}
		else {
			for (int t=0; t<inst->T; t++, c++) {	// Note: c is iterated in the secondary-loops
				pi_T[g][t] += inst->sceProb[s] * duals[c] * ( inst->capacity[g] );
			}
		}
	}
	
	// production amounts
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++, c++) {
			pi_T[g][t] += inst->sceProb[s] * duals[c] * inst->min_prod_lim[g];
		}
	}
	
	// ramp up constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=1; t<inst->T; t++, c++) {
			pi_b += inst->sceProb[s] * duals[c] * inst->ramp_u_lim[g];
		}
	}
	
	// ramp down constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=1; t<inst->T; t++, c++) {
			pi_b += inst->sceProb[s] * duals[c] * inst->ramp_d_lim[g];
		}
	}
	
	// model-specific constraints
	if (model_id == 0)
	{
		// aggregated demand constraints
		for (int t=0; t<inst->T; t++, c++) {
			pi_b += inst->sceProb[s] * duals[c] * inst->aggDemand[t];
		}
	}
	else if (model_id == 1)
	{
		// Flow upper bounds
		for (int l=0; l<inst->L; l++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += inst->sceProb[s] * duals[c] * inst->max_arc_flow[l];
			}
		}

		// Flow lower bounds
		for (int l=0; l<inst->L; l++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += inst->sceProb[s] * duals[c] * inst->min_arc_flow[l];
			}
		}
		
		// Phase-angle upper bounds
		for (int b=0; b<inst->B; b++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += inst->sceProb[s] * duals[c] * (2*inst->pi);
			}
		}

		// DC-approximation to AC flow
		// - skipped, as all necessary coefs are 0
		// - important: these constraints are not added to cons array

		// Flow balance
		for (int b=0; b<inst->B; b++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += inst->sceProb[s] * duals[c] * inst->demand[b][t];
			}
		}

	}
	else {
		cout << "!! Error: Unknown model id" << endl;
	}
}

void subprob::reset_cut_coefs () {
	pi_b = 0.0;
	for (int g=0; g<inst->G; g++) {
		std::fill(pi_T[g].begin(), pi_T[g].end(), 0.0);
	}
}

double subprob::getRecourseObjValue() {
	return recourse_obj_val;
}

void subprob::get_feasibility_cut_coefs(int &s)
{
	// reset earlier (optimality-)cut coefficient entries
	reset_cut_coefs();
	
	// get the extreme ray (based on the algorithm used to solve the subproblem)
	if ( cplex.getParam(IloCplex::RootAlg) == IloCplex::Dual ) {
	
		// Note: dualFarkas function mixes up the order of constraints
		//       We initialize an empty farkas constraint set, an empty dual-multiplier
		//		 vector, and a map to retrieve the correct values.
		IloRangeArray	farkasCons (env);
		
		// get the ray
		cplex.dualFarkas(farkasCons, duals);		// (duals is the ray...)
		
		// initialize the map
		farkasMap.clear();
		for (int c=0; c<farkasCons.getSize(); c++) {
			farkasMap.insert( pair<IloInt, double> (farkasCons[c].getId(), duals[c]) );
		}
	}
	else {
		cout << "!!Error: Cannot retrieve the extreme ray" << endl;
		exit(1);
	}
	
	// iterate over the constraints to write the Benders' feasibility cut
	int c=0;
	
	// capacity constraints
	for (int g=0; g<inst->G; g++) {
		map<int, vector<vector<double> > >::iterator it = inst->rndSupply.find(g);
		
		if ( it != inst->rndSupply.end() ) {
			for (int t=0; t<inst->T; t++, c++) {	// Note: c is iterated in the secondary-loops
				pi_T[g][t] += farkasMap[ cons[c].getId() ] * ( inst->rndSupply[g][s][t] );
			}
		}
		else {
			for (int t=0; t<inst->T; t++, c++) {	// Note: c is iterated in the secondary-loops
				pi_T[g][t] += farkasMap[ cons[c].getId() ] * ( inst->capacity[g] );
			}
		}
	}
	
	// production amounts
	for (int g=0; g<inst->G; g++) {
		for (int t=0; t<inst->T; t++, c++) {
			pi_T[g][t] += farkasMap[ cons[c].getId() ] * inst->min_prod_lim[g];
		}
	}
	
	// ramp up constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=1; t<inst->T; t++, c++) {
			pi_b += farkasMap[ cons[c].getId() ] * inst->ramp_u_lim[g];
		}
	}
	
	// ramp down constraints
	for (int g=0; g<inst->G; g++) {
		for (int t=1; t<inst->T; t++, c++) {
			pi_b += farkasMap[ cons[c].getId() ] * inst->ramp_d_lim[g];
		}
	}
	
	// model-specific constraints
	if (model_id == 0)
	{
		// aggregated demand constraints
		for (int t=0; t<inst->T; t++, c++) {
			pi_b += farkasMap[ cons[c].getId() ] * inst->aggDemand[t];
		}
	}
	else if (model_id == 1)
	{
		// Flow upper bounds
		for (int l=0; l<inst->L; l++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += farkasMap[ cons[c].getId() ] * inst->max_arc_flow[l];
			}
		}
		
		// Flow lower bounds
		for (int l=0; l<inst->L; l++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += farkasMap[ cons[c].getId() ] * inst->min_arc_flow[l];
			}
		}
		
		// Phase-angle upper bounds
		for (int b=0; b<inst->B; b++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += farkasMap[ cons[c].getId() ] * (2*inst->pi);
			}
		}
		
		// DC-approximation to AC flow
		// - skipped, as all necessary coefs are 0
		// - important: these constraints are not added to cons array
		
		// Flow balance
		for (int b=0; b<inst->B; b++) {
			for (int t=0; t<inst->T; t++, c++) {
				pi_b += farkasMap[ cons[c].getId() ] * inst->demand[b][t];
			}
		}
		
	}
	else {
		cout << "!! Error: Unknown model id" << endl;
	}
}

double subprob::computeLowerBound ()
{
    double bound = 0;

    IloEnv          env;
    IloModel        model   (env);
    IloRangeArray   cons    (env);
    IloCplex        cplex   (env);
    IloExpr         obj_func(env);

    /** Formulate the bounding problem **/

    /** Decision Variables **
    IloArray< IloNumVarArray > x (env, inst->G);    // generator states
    IloArray< IloNumVarArray > s (env, inst->G);    // start up
    IloArray< IloNumVarArray > z (env, inst->G);    // shut down
    IloArray< IloNumVarArray > p (env, inst->G);	// production amounts
    for (int g=0; g<inst->G; g++) {
        x[g] = IloNumVarArray(env, inst->T, 0, 1, ILOFLOAT);
        s[g] = IloNumVarArray(env, inst->T, 0, 1, ILOFLOAT);
        z[g] = IloNumVarArray(env, inst->T, 0, 1, ILOFLOAT);
        p[g] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
    }
    
    /** Constraints **
    // state constraints
    for (int g=0; g<inst->G; g++) {
        
        // t=0: generators are assumed to be turned on
        model.add( x[g][0] - 1 == s[g][0] - z[g][0] );
        
        // t>0
        for (int t=1; t<inst->T; t++) {
            model.add( x[g][t] - x[g][t-1] == s[g][t] - z[g][t] );
        }
    }
    
    // minimum uptime/downtime constraints
    for (int g=0; g<inst->G; g++)
    {
        // turn on inequalities
        for (int t=1; t<=inst->T; t++)
        {
            IloExpr lhs (env);
            for (int i = t-inst->min_u_time[g]+1; i<=t; i++)
            {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += 0;	// otherwise, generator is assumed to be operational (but turned on way earlier in the past)
            }
            model.add( lhs <= x[g][t-1] );
            lhs.end();
        }
    }
    
    for (int g=0; g<inst->G; g++)
    {
        // turn off inequalities
        for (int t=1; t<=inst->T; t++)
        {
            IloExpr lhs (env);
            for (int i = t-inst->min_d_time[g]+1; i<=t; i++)  {
                if (i-1 >= 0)	lhs += s[g][i-1];
                else			lhs += 0;	// otherwise, generator is assumed to be operational (but turned on way earlier in the past)
            }
            
            if (t-inst->min_d_time[g]-1 >= 0)	model.add( lhs <= 1 - x[g][t-inst->min_d_time[g]-1] );
            else								model.add( lhs <= 1 - 1 );	// assumed to be remaining on for a long, long, while
            lhs.end();
        }
    }

    // capacity constraints
    for (int g=0; g<inst->G; g++) {
        for (int t=0; t<inst->T; t++) {
            IloRange con;
            if (inst->must_run[g])  { con = IloRange (env, 0, p[g][t] - x[g][t] * inst->capacity[g], 0); }
            else                    { con = IloRange (env, -IloInfinity, p[g][t] - x[g][t] * inst->capacity[g], 0); }
            cons.add( con );
        }
    }

    // minimum production constraints
    for (int g=0; g<inst->G; g++) {
        for (int t=0; t<inst->T; t++) {
            model.add( p[g][t] >= inst->min_prod_lim[g] * x[g][t] );
        }
    }
    
    // ramping constraints
    for (int g=0; g<inst->G; g++) {
        for (int t=1; t<inst->T; t++) {
            model.add( p[g][t] - p[g][t-1] <= inst->ramp_u_lim[g] );
            model.add( p[g][t-1] - p[g][t] <= inst->ramp_d_lim[g] );
        }
    }
    
    // model-specific constraints, and the objective functions
    if (model_id == 0)
    {
        IloNumVarArray L (env, inst->T, 0, IloInfinity, ILOFLOAT);	// load shedding
        
        // aggregated demand constraints
        for (int t=0; t<inst->T; t++) {
            IloExpr expr (env);
            for (int g=0; g<inst->G; g++)	expr += p[g][t];
            expr += L[t];
            model.add( expr >= inst->aggDemand[t] );
        }
        
        // objective function
        for (int g=0; g<inst->G; g++) {
            for (int t=0; t<inst->T; t++) {
                obj_func += inst->var_cost[g] * p[g][t];
            }
        }
        
        for (int t=0; t<inst->T; t++) {
            obj_func += inst->load_shedding_penalty * L[t];
        }
    }
    else if (model_id == 1)
    {
        // node-based formulation's decision vars
        IloArray< IloNumVarArray > L (env, inst->B);        // load-shedding
        IloArray< IloNumVarArray > T (env, inst->B);		// phase angles
        for (int b=0; b<inst->B; b++) {
            L[b] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
            T[b] = IloNumVarArray(env, inst->T, 0, IloInfinity, ILOFLOAT);
        }
        
        IloArray< IloNumVarArray > F (env, inst->L);		// flows
        for (int l=0; l<inst->L; l++) {
            F[l] = IloNumVarArray(env, inst->T, -IloInfinity, IloInfinity, ILOFLOAT);
        }
        
        // Flow upper and lower bounds
        for (int l=0; l<inst->L; l++) {
            for (int t=0; t<inst->T; t++) {
                model.add( inst->min_arc_flow[l] <= F[l][t] <= inst->max_arc_flow[l] );
            }
        }
        
        // Phase-angle upper bounds
        for (int b=0; b<inst->B; b++) {
            for (int t=0; t<inst->T; t++) {
                model.add( T[b][t] <= 2*inst->pi );
            }
        }
        
        // DC-approximation to AC flow
        for (int l=0; l<inst->L; l++) {
            for (int t=0; t<inst->T; t++) {
                int orig = inst->arcs[l].first;
                int dest = inst->arcs[l].second;
                
                model.add( F[l][t] == inst->susceptance[l] * (T[orig][t] - T[dest][t]) );
            }
        }
        
        // Flow balance
        for (int b=0; b<inst->B; b++) {
            for (int t=0; t<inst->T; t++) {
                
                IloExpr expr (env);
                
                // production
                for (int g=0; g<inst->G; g++) {
                    if (inst->generator_loc[g] == b) expr += p[g][t];
                }
                
                // incoming arcs
                for (int l=0; l<inst->L; l++) {
                    if (inst->arcs[l].second == b) expr += F[l][t];
                }
                
                // outgoing arcs
                for (int l=0; l<inst->L; l++) {
                    if (inst->arcs[l].first == b) expr -= F[l][t];
                }
                
                // load shedding
                expr += L[b][t];
                
                // constraint
                model.add( expr == inst->demand[b][t] );

                // free up memory
                expr.end();
            }
        }
        
        // objective function
        for (int g=0; g<inst->G; g++) {
            for (int t=0; t<inst->T; t++) {
                obj_func += inst->var_cost[g] * p[g][t];
            }
        }
        
        for (int b=0; b<inst->B; b++) {
            for (int t=0; t<inst->T; t++) {
                obj_func += inst->load_shedding_penalty * L[b][t];
            }
        }
    }
    else {
        cout << "!! Error: Unknown model id" << endl;
    }

    model.add( IloMinimize(env, obj_func) );
    model.add( cons );
    
    cplex.extract(model);
	/* TODO: Commented by HG after encountering "Invalid arguments candidates are: void setOut(? &)
    cplex.setOut(env.getNullStream()); */
    
    /** Change model components according to scenario, and solve the subproblems **
    for (int s=0; s<inst->S; s++) {
        
        // setup subproblem s
        for (int g=0, c=0; g<inst->G; g++) {    // only the capacity constraints contain randomness
            
            map<int, vector<vector<double> > >::iterator it = inst->rndSupply.find(g);
            
            if ( it != inst->rndSupply.end() ) {
                for (int t=0; t<inst->T; t++, c++) {	// Note: c is iterated in the secondary-loop
                    cons[c].setLinearCoef(x[g][t], -inst->rndSupply[g][s][t]);
                }
            }
            else {
                c += inst->T;
            }
        }    

        // solve the subproblem
        bool status = cplex.solve();
        if (!status) {
            cout << "Subproblem " << s << " returned the status " << cplex.getCplexStatus() << endl;
            exit(1);
        }
        
        bound += inst->sceProb[s] * cplex.getObjValue();
    }
    env.end();
    
    cout << "Computed recourse lower bound = " << bound << endl;
    
    return bound;
}
	 */
