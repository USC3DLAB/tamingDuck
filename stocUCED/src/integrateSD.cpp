/*
 * integrateSD.c
 *
 *  Created on: Dec 22, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *
 * Please send your comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "integrateSD.hpp"

long long	MEM_USED = 0;	/* Amount of memory allocated each iteration */
stringC	outputDir;			/* output directory */
configType	config;			/* algorithm tuning parameters */

int integrateSD(instance &inst, EDmodel &rtED, string probName, string &configPath, ScenarioType stocObserv, int t0, double &objVal) {
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;
	vectorC edSols;

	int j=0;
	
	/* TODO: Create a common output directory: Setup outputDir */
	outputDir = (stringC) malloc(BLOCKSIZE*sizeof(char));
	strcpy(outputDir, "./");

	/* read algorithm configuration files */
	if ( readConfig(configPath) ) {
		perror("failed to read algorithm configuration file");
		goto TERMINATE;
	}

	/* Build the oneProblem structure */
	if ( (orig = buildOneProblem(rtED.model, rtED.cplex, "rted")) == NULL ) {
		perror("Failed to build the oneProblem structure for 2SD.\n");
		goto TERMINATE;
	}

#if defined(WRITE_PROB)
	testIntegrateSD(orig);
#endif

	/* Build the timeType structure */
	if ( (tim = buildTimeType(rtED.model, rtED.timeRows, rtED.timeCols)) == NULL) {
		perror("Failed to build the timeType structure for 2SD.\n");
		goto TERMINATE;
	}

	/* Build the stocType structure */
	if ( (stoc = buildStocType(rtED.model, rtED.stocRows, stocObserv, t0)) == NULL) {
		perror("Failed to build the stocType structure for 2SD.\n");
		goto TERMINATE;
	}

	/* launch the algorithm */
	if ( algo(orig, tim, stoc, (stringC) probName.c_str(), &edSols, &objVal) ) {
		perror("failed to solve the problem using 2SD");
		goto TERMINATE;
	}

	/* extract the solution */
	// Important: Variable declaration order must not be altered.
	j=0;
	for (int g=0; j<inst.powSys->numGen*2; g++, j=j+2) {
		inst.solution.usedGen_ED[g][t0] = edSols[j];
		inst.solution.overGen_ED[g][t0] = edSols[j+1];
		inst.solution.g_ED[g][t0]		= edSols[j] + edSols[j+1];
	}
	
	for (int b=0; j<inst.powSys->numGen*2+inst.powSys->numBus*3; b++, j=j+3) {
		inst.solution.loadShed_ED[b][t0] = edSols[j+1];
	}

	/* free memory */
	mem_free(edSols);
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);
	closeSolver();
	return 0;

	TERMINATE:
	mem_free(edSols);
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);
	closeSolver();
	return 1;
}//END integrateSD()

int readConfig(string &configPath) {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status;

	fptr = fopen( (configPath + "config.sd").c_str() , "r");
	if ( fptr == NULL ) {
		perror("failed to open configuration file");
		return 1;
	}

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "RUN_SEED")))
			fscanf(fptr, "%lld", &config.RUN_SEED);
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MASTERTYPE")))
			fscanf(fptr, "%d", &config.MASTERTYPE);
		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "TAU")))
			fscanf(fptr, "%d", &config.TAU);
		else if (!(strcmp(line, "MIN_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MIN_QUAD_SCALAR);
		else if (!(strcmp(line, "MAX_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MAX_QUAD_SCALAR);
		else if (!(strcmp(line, "R1")))
			fscanf(fptr, "%lf", &config.R1);
		else if (!(strcmp(line, "R2")))
			fscanf(fptr, "%lf", &config.R2);
		else if (!(strcmp(line, "R3")))
			fscanf(fptr, "%lf", &config.R3);
		else if (!(strcmp(line, "PI_EVAL_START")))
			fscanf(fptr, "%d", &config.PI_EVAL_START);
		else if (!(strcmp(line, "PI_CYCLE")))
			fscanf(fptr, "%d", &config.PI_CYCLE);
		else if (!(strcmp(line, "SCAN_LEN")))
			fscanf(fptr, "%d", &config.SCAN_LEN);
		else if (!(strcmp(line, "EVAL_FLAG")))
			fscanf(fptr, "%d", &config.EVAL_FLAG);
		else if (!(strcmp(line, "EVAL_SEED")))
			fscanf(fptr, "%lld", &config.EVAL_SEED);
		else if (!(strcmp(line, "EVAL_MIN_ITER")))
			fscanf(fptr, "%d", &config.EVAL_MIN_ITER);
		else if (!(strcmp(line, "EVAL_ERROR")))
			fscanf(fptr, "%lf", &config.EVAL_ERROR);
		else if (!(strcmp(line, "PRE_EPSILON")))
			fscanf(fptr, "%lf", &config.PRE_EPSILON);
		else if (!(strcmp(line, "EPSILON")))
			fscanf(fptr, "%lf", &config.EPSILON);
		else if (!(strcmp(line, "BOOTSTRAP_REP")))
			fscanf(fptr, "%d", &config.BOOTSTRAP_REP);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			perror("unrecognized parameter in configuration file encountered while reading readConfig");
		}
	}

	fclose(fptr);

	return 0;
}//END readConfig()

oneProblem *buildOneProblem(IloModel &model, IloCplex &cplex, string probName) {
	oneProblem *orig;

	orig = (oneProblem *) malloc(sizeof(oneProblem));

	/* Set default values */
	orig->objx = orig->bdl = orig->bdu = orig->rhsx = orig->matval = NULL;
	orig->matind = orig->matbeg = orig->matcnt = NULL;
	orig->name = orig->senx = orig->ctype = orig->objname = orig->cstore = orig->rstore = NULL;
	orig->cname = orig->rname = NULL;

	orig->mac = orig->macsz = orig->mar = orig->marsz = orig->numBin = orig->numInt = orig->numnz = orig->cstorsz = orig->rstorsz = 0;

	orig->objsen = 1;

	/* Problem name and size */
	orig->name 	= (stringC) malloc(NAMESIZE*sizeof(char)); probName.copy(orig->name, probName.size(), 0);
	orig->mac   = orig->macsz = cplex.getNcols();
	orig->mar   = orig->marsz = cplex.getNrows();
	orig->numnz = orig->matsz = cplex.getNNZs();

	/* Allocate memory for column elements of SD-oneProblem */
	orig->bdl 	= (vectorC) calloc(orig->mac,sizeof(double));
	orig->bdu 	= (vectorC) calloc(orig->mac, sizeof(double));
	orig->ctype = (stringC) calloc(orig->mac, sizeof(char));
	orig->cstore = (stringC) calloc(orig->mac*NAMESIZE, sizeof(char));
	orig->cname = (stringC *) calloc(orig->mac, sizeof(stringC));

	/* Allocate memory for row elements of SD-oneProblem */
	orig->senx 	= (stringC) calloc(orig->mar, sizeof(char));
	orig->rhsx 	= (vectorC) calloc(orig->mar, sizeof(double));
	orig->rstore = (stringC) calloc(orig->mar, NAMESIZE*sizeof(char));
	orig->rname	= (stringC *) calloc(orig->mar, sizeof(stringC));

	/* Allocate memory to hold objective function coefficients. */
	IloObjective obj;
	orig->objx = (vectorC) calloc(orig->mac, sizeof(double));

	/* Prepare to read the constraint matrix */
	map<string, int> colID;
	vector<vector<int>> idx;
	vector<vector<double>> vals;
	for ( int c = 0; c < orig->mac; c++ ) {
		idx.push_back(vector<int> ());
		vals.push_back(vector<double> ());
	}
	/* Allocate memory for constraint matrix elements of SD-oneProblem */
	orig->matcnt = (intvec) malloc(orig->mac*sizeof(int));
	orig->matbeg = (intvec) malloc(orig->mac*sizeof(int));
	orig->matind = (intvec) malloc(orig->numnz*sizeof(int));
	orig->matval = (vectorC) malloc(orig->numnz*sizeof(double));

	/* Loop through the model iterates and capture objective function and variables */
	{
		IloModel::Iterator it(model);
		for (int c = 0, m = 0; it.ok(); ++it ) {
			if ( (*it).isObjective() ) {
				/* Objective function */
				obj = (*it).asObjective();
				orig->objname = (stringC) malloc(NAMESIZE*sizeof(char));
				orig->objsen = obj.getSense();
				string tempName = obj.getName();
				if(!tempName.empty()) {
					tempName = tempName + '\0';
					tempName.copy(orig->objname, tempName.size(), 0);
				}
			}
			else if ( (*it).isVariable() ) {
				/* Decision variable */
				IloNumVar var = (*it).asVariable();

				orig->bdl[c] = var.getLB();
				orig->bdu[c] = var.getUB();

				IloNumVar::Type ctype = var.getType();
				if ( ctype == 1 ) { orig->ctype[c] = 'I'; }
				else if ( ctype == 2 ) { orig->ctype[c] = 'C';}
				else if ( ctype == 3 ) { orig->ctype[c] = 'B'; }
				else { perror ("Unknown column type."); return NULL; }

				string tempName = var.getName();
				orig->cname[c] = orig->cstore + m;
				orig->cstorsz += tempName.size() + 1;
				for ( auto it2 = tempName.begin(); it2 != tempName.end()+1; ++it2 )
					orig->cstore[m++] = *it2;

				colID.insert( pair<string, int> (tempName, c) );
				c++;
			}
			else if ( (*it).isConstraint() )
				continue;
			else {
				perror("Unknown iterate type.\n");
				return NULL;
			}
		}
	}

	/* Loop through the iterate again, this time to capture the constraints. Need to loop separately to makes sure all the variables
	 * are read before the constraints are being read. This ensures that the constraint matrix is built properly. */
	{
		IloModel::Iterator it(model);
		for (int r = 0, m = 0; it.ok(); ++it) {
			if ( (*it).isConstraint() ) {
				/* Constraint */
				IloConstraint con = (*it).asConstraint();
				string tempName = con.getName();
				orig->rname[r] = orig->rstore + m;
				orig->rstorsz += tempName.size()+1;
				for ( auto it2 = tempName.begin(); it2 != tempName.end()+1; ++it2 )
					orig->rstore[m++] = *it2;

				IloRangeI *impl = dynamic_cast<IloRangeI *>(con.getImpl());
				if ( impl ) {
					IloRange rng(impl);
					for (IloExpr::LinearIterator it2 = IloExpr(rng.getExpr()).getLinearIterator(); it2.ok(); ++it2) {
						string tempName = (it2.getVar()).getName();
						double tempVal = it2.getCoef();
						auto i = colID.find(tempName);
						if ( i == colID.end() ) { perror("Could not find the column name while parsing row.\n"); return NULL; }
						idx[i->second].push_back(r);
						vals[i->second].push_back(tempVal);
					}
					IloNum lb = rng.getLB();
					IloNum ub = rng.getUB();
					if ( lb == ub ) {
						orig->senx[r] = 'E';
						orig->rhsx[r] = lb;
					}
					else if ( lb == -IloInfinity ) {
						orig->senx[r] = 'L';
						orig->rhsx[r] = ub;
					}
					else if ( ub == IloInfinity ) {
						orig->senx[r] = 'G';
						orig->rhsx[r] = lb;
					}
					else {
						cout << rng << endl;
						perror("Unknown bounds on constraints.\n");
					}
				}
				r++;
			}
		}
	}

	/* Allocate memory for objective coefficients */
	for ( IloExpr::LinearIterator it = IloExpr(obj.getExpr()).getLinearIterator(); it.ok(); ++it ) {
		auto i = colID.find((it.getVar()).getName());
		if ( i == colID.end() ) { perror("Could not find the column name while parsing objective function.\n"); return NULL; }
		orig->objx[i->second] = it.getCoef();
	}

	/* Setup constraint matrix */
	int cnt = 0;
	for ( int c = 0; c < orig->mac; c++ ) {
		orig->matcnt[c] = idx[c].size(); orig->matbeg[c] = cnt;
		for ( int n = 0; n < orig->matcnt[c]; n++ ) {
			orig->matval[cnt] = vals[c][n];
			orig->matind[cnt++] = idx[c][n];
		}
	}

	return orig;
}//END buildOneProblem()

timeType *buildTimeType(IloModel &model, vector<string> rowNames, vector<string> colNames) {
	timeType *tim = NULL;

	map<string, int> colID;
	map<string, int> rowID;

	/* Get column and row names */
	IloModel::Iterator it(model);
	for (int c = 0, r = 0; it.ok(); ++it ) {
		if ( (*it).isVariable() ) {
			/* Decision variable */
			IloNumVar var = (*it).asVariable();

			string tempName = var.getName();
			colID.insert( pair<string, int> (tempName, c) );
			c++;
		}
		else if ( (*it).isConstraint() ) {
			/* Constraints */
			IloConstraint con = (*it).asConstraint();

			string tempName = con.getName();
			rowID.insert( pair<string, int> (tempName, r) );
			r++;
		}
	}

	/* allocate memory and initialize */
	tim = (timeType *) malloc(sizeof(timeType));
	tim->probName = (stringC) malloc(NAMESIZE*sizeof(char));
	tim->stgNames = (stringC *) malloc(2*sizeof(stringC));
	tim->row = (intvec) malloc(2*sizeof(int));
	tim->col = (intvec) malloc(2*sizeof(int));
	tim->rowStg = tim->colStg = NULL;

	/* Default parameters */
	tim->type = 0;
	strcpy(tim->probName, "rted");
	tim->numStages = 2;
	tim->numCols = 0; tim->numRows = 0;

	/* Use row and column names to set the time decomposition parameters. */
	for ( int t = 0; t < (int) rowNames.size(); t++ ) {
		tim->stgNames[t] = (stringC) malloc(NAMESIZE*sizeof(char));
		sprintf(tim->stgNames[t], "stage%d", t);

		auto idx = colID.find(colNames[t]);
		if ( idx == colID.end() ) { perror("Could not find the column name while parsing the time decomposition details.\n"); return NULL; }
		tim->col[t] = idx->second;

		idx = rowID.find(rowNames[t]);
		if ( idx == rowID.end() ) { perror("Could not find the row name while parsing the time decomposition details.\n"); return NULL; }
		tim->row[t] = idx->second;
	}

	return tim;
}//END buildTimeType()

stocType *buildStocType(IloModel &model, vector<string> stocRows, ScenarioType stocObserv, int t0) {
	stocType *stoc = NULL;

	int maxOmegas = stocRows.size();

	/* allocate memory to stocType and initialize elements */
	stoc = (stocType *) calloc(1, sizeof(stocType));
	stoc->type = (stringC) malloc(NAMESIZE*sizeof(char));
	stoc->col = (intvec) calloc(maxOmegas, sizeof(int));
	stoc->row = (intvec) calloc(maxOmegas, sizeof(int));

	stoc->mean = (vectorC) calloc(maxOmegas, sizeof(double));
	stoc->numVals = (intvec) calloc(maxOmegas, sizeof(int));
	stoc->vals = (vectorC *) calloc(maxOmegas, sizeof(vectorC));
	stoc->probs = (vectorC *) calloc(maxOmegas, sizeof(vectorC));

	stoc->groupBeg = (intvec) calloc(maxOmegas, sizeof(int));
	stoc->numPerGroup = (intvec) calloc(maxOmegas, sizeof(int));

	stoc->numCipher = 0; stoc->sim = FALSE; stoc->arma = NULL;

	/* Default parameters */
	strcpy(stoc->type, "BLOCKS_DISCRETE");
	stoc->numGroups = 1; stoc->numPerGroup[0] = maxOmegas;
	stoc->numOmega = maxOmegas;

	map<string, int> colID;
	map<string, int> rowID;

	/* Get column and row names */
	IloModel::Iterator it(model);
	for (int c = 0, r = 0; it.ok(); ++it ) {
		if ( (*it).isVariable() ) {
			/* Decision variable */
			IloNumVar var = (*it).asVariable();

			string tempName = var.getName();
			colID.insert( pair<string, int> (tempName, c) );
			c++;
		}
		else if ( (*it).isConstraint() ) {
			/* Constraints */
			IloConstraint con = (*it).asConstraint();

			string tempName = con.getName();
			rowID.insert( pair<string, int> (tempName, r) );
			r++;
		}
	}

	/* Number of blocks and their probabilities are stored in index-0 */
	stoc->numVals[0] = stocObserv.cnt;
	stoc->probs[0] = (vectorC) calloc(stoc->numVals[0], sizeof(double));
	for (int cnt = 0; cnt < stoc->numVals[0]; cnt++ ) {
		stoc->probs[0][cnt] = (double) 1/ (double) stoc->numVals[0];
	}

	for ( int m = 0; m < maxOmegas; m++ ) {
		stoc->col[m] = -1;

		auto idx = rowID.find(stocRows[m]);
		if ( idx == rowID.end() ) { perror("Could not find the row name while parsing the time decomposition details.\n"); return NULL; }
		stoc->row[m] = idx->second;

		stoc->mean[m] = 0;

		stoc->vals[m] = (vectorC) calloc(stoc->numVals[0], sizeof(double));
		for ( int cnt = 0; cnt < stoc->numVals[0]; cnt++ ) {
			int t = m/stocObserv.numVars;
			int n = m % stocObserv.numVars;
			stoc->vals[m][cnt] = stocObserv.vals[cnt][t0+1+t][n]; //TODO: Verify
		}
	}

	return stoc;
}//END buildStocType()

oneProblem *buildOneProblem_file(string probName) {
	oneProblem *orig;
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(env);
	IloObjective   obj;
	IloNumVarArray var(env);
	IloRangeArray  rng(env);

	cplex.importModel(model, probName.c_str(), obj, var, rng);
	cplex.extract(model);

	orig = (oneProblem *) malloc(sizeof(oneProblem));

	/* Set default values */
	orig->objx = orig->bdl = orig->bdu = orig->rhsx = orig->matval = NULL;
	orig->matind = orig->matbeg = orig->matcnt = NULL;
	orig->name = orig->senx = orig->ctype = orig->objname = orig->cstore = orig->rstore = NULL;
	orig->cname = orig->rname = NULL;

	orig->mac = orig->macsz = orig->mar = orig->marsz = orig->numBin = orig->numInt = orig->numnz = orig->cstorsz = orig->rstorsz = 0;

	orig->objsen = 1; orig->type = PROB_LP;

	/* Problem name and size */
	orig->name 	= (stringC) malloc(NAMESIZE*sizeof(char)); probName.copy(orig->name, probName.size(), 0);
	orig->mac   = orig->macsz = cplex.getNcols();
	orig->mar   = orig->marsz = cplex.getNrows();
	orig->numnz = cplex.getNNZs();

	/* Objective function sense and coefficients */
	orig->objname = (stringC) malloc(NAMESIZE*sizeof(char));
	orig->objsen = obj.getSense();
	{ string tempName = obj.getName(); if(!tempName.empty()) {	tempName = tempName + '\0'; tempName.copy(orig->objname, tempName.size(), 0);}}

	/* Column information */
	int c = 0, m = 0;
	map<string, int> colID;
	/* Allocate memory for column elements of SD-oneProblem */
	orig->objx 	= (vectorC) malloc(orig->mac*sizeof(double));
	orig->bdl 	= (vectorC) malloc(orig->mac*sizeof(double));
	orig->bdu 	= (vectorC) malloc(orig->mac*sizeof(double));
	orig->ctype = (stringC) malloc(orig->mac*sizeof(char));
	orig->cstore = (stringC) malloc(orig->mac*NAMESIZE*sizeof(char));
	orig->cname = (stringC *) malloc(orig->mac*sizeof(stringC));

	for ( IloExpr::LinearIterator it = IloExpr(obj.getExpr()).getLinearIterator(); it.ok(); ++it ) {
		orig->objx[c] = it.getCoef();
		orig->bdl[c] = var[c].getLB();
		orig->bdu[c] = var[c].getUB();
		orig->ctype[c] = var[c].getType();

		string tempName = var[c].getName();
		orig->cname[c] = orig->cstore + m;
		orig->cstorsz += tempName.size()+1;
		for ( auto it = tempName.begin(); it != tempName.end()+1; ++it )
			orig->cstore[m++] = *it;

		colID.insert( pair<string, int> (tempName, c) );
		c++;
	}

	/* Allocate memory for row elements of SD-oneProblem */
	orig->senx 	= (stringC) malloc(orig->mar*sizeof(char));
	orig->rhsx 	= (vectorC) malloc(orig->mar*sizeof(double));
	orig->rstore = (stringC) malloc(orig->mar*NAMESIZE*sizeof(char));
	orig->rname	= (stringC *) malloc(orig->mar*sizeof(stringC));

	/* Prepare to read the constraint matrix */
	vector<vector<int>> idx;
	vector<vector<double>> vals;
	for ( int c = 0; c < orig->mac; c++ ) {
		idx.push_back(vector<int> ());
		vals.push_back(vector<double> ());
	}

	/* Row information */
	m = 0;
	for ( int r = 0; r < orig->mar; r++ ) {
		for (IloExpr::LinearIterator it2 = IloExpr(rng[r].getExpr()).getLinearIterator(); it2.ok();++it2) {
			string tempName = (it2.getVar()).getName();
			orig->rname[r] = orig->rstore + m;
			orig->rstorsz += tempName.size()+1;
			for ( auto it2 = tempName.begin(); it2 != tempName.end()+1; ++it2 )
				orig->rstore[m++] = *it2;

			double tempVal = it2.getCoef();
			auto i = colID.find(tempName);
			if ( i == colID.end() ) { perror("Could not find the column name while parsing row.\n"); return NULL; }
			idx[i->second].push_back(r);
			vals[i->second].push_back(tempVal);
		}

		IloNum lb = rng[r].getLB();
		IloNum ub = rng[r].getUB();
		if ( lb == ub ) {
			orig->senx[r] = 'E';
			orig->rhsx[r] = lb;
		}
		else if ( lb == -IloInfinity ) {
			orig->senx[r] = 'L';
			orig->rhsx[r] = ub;
		}
		else if ( ub == IloInfinity ) {
			orig->senx[r] = 'G';
			orig->rhsx[r] = lb;
		}
	}

	/* Allocate memory for constraint matrix elements of SD-oneProblem */
	orig->matcnt = (intvec) malloc(orig->mac*sizeof(int));
	orig->matbeg = (intvec) malloc(orig->mac*sizeof(int));
	orig->matind = (intvec) malloc(orig->numnz*sizeof(int));
	orig->matval = (vectorC) malloc(orig->numnz*sizeof(double));

	/* Setup constraint matrix */
	int cnt = 0;
	for ( int c = 0; c < orig->mac; c++ ) {
		orig->matcnt[c] = idx[c].size(); orig->matbeg[c] = cnt;
		for ( int n = 0; n < orig->matcnt[c]; n++ ) {
			orig->matval[cnt] = vals[c][n];
			orig->matind[cnt++] = idx[c][n];
		}
	}

	return orig;
}//END buildOneProblem()

void testIntegrateSD(oneProblem *orig) {
	LPptr lp;

	openSolver();

	char fname[BLOCKSIZE];
	if ( !(lp = setupProblem(orig->name, orig->type, orig->mac, orig->mar, orig->objsen, orig->objx, orig->rhsx, orig->senx,
			orig->matbeg, orig->matcnt, orig->matind, orig->matval, orig->bdl, orig->bdu, NULL, orig->cname, orig->rname, orig->ctype)) ) {
		perror("failed to setup stage problem in solver");
	}
	sprintf(fname, "builtOneProblem.lp");
	if ( writeProblem(lp, fname) ) {
		perror("failed to write stage problem");
	}

	closeSolver();

}//END testIntegrateSD()
