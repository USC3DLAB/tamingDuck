/*
 * solver.c
 *
 *  Created on: Apr 20, 2014
 *      Author: Harsha Gangammanavar
 */

#include "utils.h"
#include "solver.h"

extern stringC outputDir;
ENVptr	env;

int solveProblem(LPptr lp, stringC pname, int type, int *status) {
	int		aggres = 0;

	solveagain:
	switch  ( type ) {
	case PROB_LP:
        changeSolverType(ALG_AUTOMATIC);
		(*status) = CPXlpopt(env, lp);
		break;
	case PROB_QP:
        changeSolverType(ALG_CONCURRENT);
		(*status) = CPXbaropt(env, lp);
		break;
	case PROB_MILP:
		(*status) = CPXmipopt(env, lp);
		break;
	case PROB_MIQP:
		(*status) = CPXmipopt(env, lp);
		break;
	default:
		break;
	}

RESOLVE:
	(*status) = CPXgetstat(env, lp);
	if ((*status) != STAT_OPTIMAL && (*status) != MIP_OPTIMAL && (*status)!= MIP_OPTIMAL_TOL  ) {
		if ((*status) == STAT_INFEASIBLE || (*status) == MIP_INFEASIBLE) {
			writeProblem(lp, "error.lp");
			return (*status);
		}
		else if ( (type == PROB_LP || type == PROB_QP) && (*status) == 6 ){
			setIntParam(PARAM_SCAIND, 1);
			aggres++;
			if(aggres == 1)
				goto solveagain;
			else
				goto skip;
        } else if ( (type == PROB_QP) && (*status) == 2 ) {
            changeSolverType(ALG_PRIMAL);
            (*status) = CPXlpopt(env, lp);
            goto RESOLVE;
        }
        else {
			solverErrmsg((*status));
			writeProblem(lp, "error.lp");
			return (*status);
		}
	}

	skip: if (aggres)
		setIntParam(PARAM_SCAIND, 0);

	return 0;
}//END solveProblem()

int getProbType(LPptr lp) {
	int status;

	status = CPXgetprobtype (env, lp);
	if ( status < -1 ) {
		solverErrmsg(status);
		exit(1);
	}
	else
		return status;
}//END

double getObjective(LPptr lp, int type) {
	double ans;
	int status;

	/* Get the objective function value only */
	CPXsolution(env, lp, &status, &ans, NULL, NULL, NULL, NULL);
	if ( type == PROB_LP || type == PROB_QP ) {
		if(status != STAT_OPTIMAL && status != 6) {
			solverErrmsg(status);
			exit(1);
		}
	}
	else
		if(status != MIP_OPTIMAL && status != MIP_OPTIMAL_TOL ) {
			solverErrmsg(status);
			exit(1);
		}

	return ans;
}//END getObjective()

int getPrimal(LPptr lp, vectorC X, int length) {
	int status;
    
	status = CPXgetx(env, lp, X+1, 0, length-1);
	if ( status )
		solverErrmsg(status);
	else
		X[0] = oneNorm(X+1, length);

	return status;
}//END getPrimal()

double getPrimalPoint(LPptr lp, int idx) {
	int 	status;
	double 	x[1];

	status = CPXgetx(env, lp, x, idx, idx);
	if ( status )
		solverErrmsg(status);

	return x[0];
}//END getPrimalPoint

int getDual(LPptr lp, vectorC Pi, int length) {
	int	status;

	status = CPXgetpi(env, lp, Pi+1, 0, length-1);
	if(status)
		solverErrmsg(status);
	else
		Pi[0] = oneNorm(Pi+1, length);

	return status;
}//END getDual()

int getDualSlacks(LPptr lp, vectorC Dj, int length){
	int	status;

	status = CPXgetdj(env, lp, Dj+1, 0, length-1);
	if(status)
		solverErrmsg(status);
	else
		Dj[0] = oneNorm(Dj+1, length);

	return status;
}//END getDualSlacks()

int getBasis(LPptr lp, intvec cstat, intvec rstat){
	int status;

	status = CPXgetbase(env, lp, cstat, rstat);
	if ( status )
		solverErrmsg(status);

	return status;
}//END get_basis()

void openSolver(){
	int 	status = 0;

	env = CPXopenCPLEX(&status);
	if ( env == NULL ) {
		solverErrmsg(status);
		errMsg("solver", "openSolver", "could not open CPLEX environment", 0);
		goto TERMINATE;
	}

	status = setIntParam(PARAM_SCRIND, OFF);
	if ( status ) {
		errMsg("solver", "open_solver", "screen output", 0);
		goto TERMINATE;
	}

	TERMINATE:
	if ( status )
		closeSolver();

}// openSolver()

void closeSolver(){
	int status;

	if (env != NULL) {
		status = CPXcloseCPLEX(&env);
		if ( status ) {
			solverErrmsg(status);
			errMsg("solver", "close_solver", "could not close CPLEX environment", 1);
		}
	}

}//END closeSolver()

int changeSolverType(int method) {
	int status = 0;

	status = setIntParam(PARAM_LPMETHOD, method);
	if (status) {
		solverErrmsg(status);
		return 1;
	}
	return 0;
}

int setIntParam(int paramname, int paramvalue){

	return CPXsetintparam (env, paramname, paramvalue);

}//END setIntParam()

void solverErrmsg(int status){
	char  errmsg[CPXMESSAGEBUFSIZE];

	CPXgeterrorstring (env, status, errmsg);
	fprintf (stderr, "%s", errmsg);

}//END solverErrmsg()

int createProblem(char *probname, LPptr *lp) {
	int status;

	(*lp) = CPXcreateprob (env, &status, probname);
	if ( (*lp) == NULL )
		solverErrmsg(status);

	return status;
}// END createProb()

int readProblem(stringC probpath, LPptr lp) {
	int status;

	status = CPXreadcopyprob (env, lp, probpath, "MPS");
	if( status ) {
		status = CPXreadcopyprob(env, lp, probpath,"LP");
		if ( status )
			solverErrmsg(status);
	}

	return status;
}// END readProb()

LPptr setupProblem(stringC name, int type, int numcols, int numrows, int objsense, vectorC objx, vectorC rhsx, stringC sense, intvec matbeg, intvec matcnt,
		intvec matind, vectorC matval, vectorC lb, vectorC ub, vectorC rngval, stringC *colname, stringC *rowname, stringC ctype) {
	LPptr	lp = NULL;
	int		*indices, c;

	createProblem(name, &lp);
	if ( lp == NULL ) {
		errMsg("solver", "setupProblem", name, 1);
		return NULL;
	}

	if ( loadProbwNames(lp, numcols, numrows, objsense, objx, rhsx, sense, matbeg, matcnt, matind, matval, lb, ub, rngval, colname, rowname) ) {
		errMsg("solver", "setupProblem", "failed to load problem onto the solver", 0);
		return NULL;
	}

	if ( type == PROB_MILP || type == PROB_MIQP ) {
		if ( !(indices = (intvec) arr_alloc(numcols, int)) )
			errMsg("allocation", "setupProblem", "indices", 0);
		for ( c = 0; c < numcols; c++ )
			indices[c] = c;
		if ( changeCtype(lp, numcols, indices, ctype) ) {
			errMsg("solver", "setupProblem", "failed to change column type", 0);
			return NULL;
		}
		mem_free(indices);
	}

	return lp;
}//END setupProblem()

int loadProblem(LPptr lp, int numcols, int numrows, int objsense, vectorC objx, vectorC rhsx, stringC sense, intvec matbeg, intvec matcnt,
		intvec matind, vectorC matval, vectorC lb, vectorC ub, vectorC rngval) {
	int status;
	status = CPXcopylp(env, lp, numcols, numrows, objsense, objx, rhsx, sense, matbeg, matcnt, matind, matval, lb, ub, rngval);
	if(status)
		solverErrmsg(status);

	return status;
}// END loadwnames_prob()

int loadProbwNames(LPptr lp, int numcols, int numrows, int objsense, vectorC objx, vectorC rhsx, stringC sense, intvec matbeg, intvec matcnt,
		intvec matind, vectorC matval, vectorC lb, vectorC ub, vectorC rngval, stringC *colname, stringC *rowname) {
	int status;
	status = CPXcopylpwnames (env, lp, numcols, numrows, objsense, objx, rhsx, sense, matbeg, matcnt, matind, matval, lb, ub, rngval, colname, rowname);
	if(status)
		solverErrmsg(status);

	return status;
}// END loadwnames_prob()

int writeProblem(LPptr lp, stringC filename) {
	int status;
	char buffer[2*BLOCKSIZE];

	strcpy(buffer,outputDir);
	strcat(buffer,filename);

	status = CPXwriteprob (env, lp, buffer, "LP");
	if (status)
		solverErrmsg(status);

	return status;
}// END write_prob()

LPptr cloneProblem(LPptr origLp) {
	int status = 0;
	LPptr cloneLp = NULL;

	cloneLp = CPXcloneprob(env, origLp, &status);
	if (status != 0 || cloneLp == NULL)
		solverErrmsg(status);

	return cloneLp;
}//END cloneProb()

int getProbName(LPptr lp, stringC probName, int len) {
	int status, surplus;

	status = CPXgetprobname (env, lp, probName, len, &surplus);
	if ( surplus < 0)
		solverErrmsg(status);

	return status;
}//END getProbName()

int getObjSen(LPptr lp) {
	int status;

	status = CPXgetobjsen(env, lp);
	if ( !(status) )
		solverErrmsg(status);

	return status;
}// END getObjsen()

int getNumRows(LPptr lp) {
	int status;

	status = CPXgetnumrows(env,lp);
	if ( !(status) )
		solverErrmsg(status);

	return status;
}// END getNumRows()

int getNumCols(LPptr lp) {
	int status;

	status = CPXgetnumcols(env, lp);
	if ( !(status) )
		solverErrmsg(status);

	return status;
}// END getNumCols()

int getNumBinary(LPptr lp){

	return CPXgetnumbin(env, lp);

}//END getNumBinary()

int getNumInt(LPptr lp){

	return CPXgetnumint(env, lp);

}//END getNumBinary()

int getCtype(LPptr lp, int start, int end, stringC ctype) {
	int status;

	status = CPXgetctype (env, lp, ctype, start, end-1);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getCtype()

int getNumnz(LPptr lp) {

	return CPXgetnumnz (env, lp);

}//END getNumnz

int getObjx(LPptr lp, int start, int end, vectorC obj) {
	int status;

	if( (status = CPXgetobj (env, lp, obj, start, end-1)) )
		solverErrmsg(status);

	return status;
}// END getObjx()

int getRhsx(LPptr lp, int start, int end, vectorC rhs) {
	int status;

	if( (status = CPXgetrhs (env, lp, rhs, start, end-1)))
		solverErrmsg(status);

	return status;
}// END getRhsx()

int getSense(LPptr lp, int start, int end, stringC sense) {
	int status;

	if ( (status = CPXgetsense (env, lp, sense, start, end-1)))
		solverErrmsg(status);

	return status;
}// END getSense()

int getCols(LPptr lp, int start, int end, intvec cmatbeg, intvec cmatind, vectorC cmatval, int cmatspace) {
	int status, nzcnt, surplus;

	status = CPXgetcols (env, lp, &nzcnt, cmatbeg, cmatind, cmatval, cmatspace, &surplus, start, end-1);
	if ( status != 0 || surplus < 0 ) {
		fprintf(stderr, "Surplus = %d", surplus);
		return 1;
	}
	if (nzcnt > cmatspace)
		printf("The allocated size and used size for constraint elements do not match.\n");

	return 0;
}// END getCols()

/* function reads the lower bound on variables indexed from _start_ to _end_ of problem input and stores them in _lb_. Returns 0 on success and 1 on
 * failure. */
int getLb(LPptr lp, int start, int end, vectorC lb) {
	int status;

	if ( (status = CPXgetlb(env, lp, lb, start, end - 1)) )
		solverErrmsg(status);

	return status;
}//END getLb()

/* function reads the upper bound on variables indexed from _start_ to _end_ of problem input and stores them in _lb_. Returns 0 on success and 1 on
 * failure. */
int getUb(LPptr lp, int start, int end, vectorC ub) {
	int status;

	if ( (status = CPXgetub(env, lp, ub, start, end - 1)))
		solverErrmsg(status);

	return status;
}//END getUb()

/* read problem objective function name into the string _objname. Returns 0 on success and 1 on failure.*/
int getObjName(LPptr lp, stringC objname) {

	return CPXgetobjname (env, lp, objname, NAMESIZE, NULL);

}// END getObjName()

int getCstoreSize(LPptr lp, int start, int end) {
	int status, surplus;

	status = CPXgetcolname (env, lp, NULL, NULL, 0, &surplus, 0, end - 1);
	if (( status != CPXERR_NEGATIVE_SURPLUS ) && ( status != 0 ))
		return 1;

	return surplus;
}// END getCstoreSize()

int getColName(LPptr lp, int start, int end, stringC *colname, stringC colnamestore, int csize) {
	int status;

	return CPXgetcolname (env, lp, colname, colnamestore, csize, &status, start, end-1);

}// END get_colname()

int getRstoreSize(LPptr lp, int start, int end) {
	int status, surplus;

	status = CPXgetrowname (env, lp, NULL, NULL, 0, &surplus, 0, end - 1);
	if (( status != CPXERR_NEGATIVE_SURPLUS ) && ( status != 0 ))
		return 1;

	return surplus;
}//END getRstoreSize()

int getRowName(LPptr lp, int start, int end, stringC *rowname, stringC rownamestore, int rsize) {
	int status;

	return CPXgetrowname (env, lp, rowname, rownamestore, rsize, &status, start, end-1);

}// END getRowName()

int getBinvC(LPptr lp, int col, vectorC a ) {
	int status;

	status = CPXbinvcol(env, lp, col, a);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBinvC

int getBhead(LPptr lp, intvec head, vectorC x) {
	int status;

	status=CPXgetbhead(env, lp, head, x);
	if ( status )
		solverErrmsg(status);

	return status;
}//END getBhead

int binvArow(LPptr lp, int i, vectorC z) {
	int status;

	status= CPXbinvarow (env, lp, i, z);
	if ( status )
		solverErrmsg(status);

	return status;
}//END binvArow

int binvAcol(LPptr lp, int i, vectorC z) {
	int status;

	status= CPXbinvacol (env, lp, i, z);
	if ( status )
		solverErrmsg(status);

	return status;
}//END binvArow

int copyQPseparable(LPptr lp, vectorC qsepvec){
	int status;
	/* NOTE: CPLEX evaluates the corresponding objective with a factor of 0.5 in front of the quadratic objective term.*/
	status = CPXcopyqpsep (env, lp, qsepvec);
	if(status)
		solverErrmsg(status);

	return status;
}//END copy_qp_separable()

int changeCoef(LPptr lp, int row, int col, double val) {
	int status;

	status = CPXchgcoef (env, lp, row, col, val);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeObjx()

int changeObjx(LPptr lp, int cnt, intvec indices, vectorC values){
	int status;

	status = CPXchgobj (env, lp, cnt, indices, values);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeObjx()

int changeRHS(LPptr lp, int cnt, intvec indices, vectorC values) {
	int status;

	status = CPXchgrhs (env, lp, cnt, indices, values);
	if(status)
		solverErrmsg(status);

	return status;
}//END changeRHS

int changeBDS(LPptr lp, int cnt, intvec indices, stringC lu, vectorC bd) {
	int status;

	status = CPXchgbds (env, lp, cnt, indices, lu, bd);
	if( status )
		solverErrmsg(status);

	return status;
}//END changeBDS

int changeCol(LPptr lp, int column, vectorC coef, int start, int stop){
	int		row, status;

	for (row = start; row < stop; row++){
		status = CPXchgcoef(env, lp, row, column, coef[row-start]);
		if( status )
			return status;
	}

	return 0;
}//END changeCol()

int changeCtype(LPptr lp, int cnt, intvec indices, stringC ctype){
	int status;

	status = CPXchgctype (env, lp, cnt, indices, ctype);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeCtype()

int changeProbType(LPptr lp, int type) {
	int status;

	status = CPXchgprobtype (env, lp, type);
	if ( status )
		solverErrmsg(status);

	return status;
}//END changeProbType()

int addRow(LPptr lp, int nzcnt, double inputRHS, char inputSense, int matbeg, intvec rmatind, vectorC rmatval) {
	static int cumul_num = 0;
	stringC	*rowname;
	char	sense[1] = {'G'};
	double	rhs[1];
	int		status, rmatbeg[1];

	rowname = arr_alloc(2, stringC);
	rowname[0] = arr_alloc(NAMESIZE, char);
	strcpy(rowname[0], "Cut    ");
	rowname[0][3] = '0' + cumul_num / 10000 % 10;
	rowname[0][4] = '0' + cumul_num / 1000 % 10;
	rowname[0][5] = '0' + cumul_num / 100 % 10;
	rowname[0][6] = '0' + cumul_num / 10 % 10;
	rowname[0][7] = '0' + cumul_num / 1 % 10;
	cumul_num++;

	rhs[0] = inputRHS;
	sense[0] = inputSense;
	rmatbeg[0] = matbeg;

	status = CPXaddrows (env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
	if ( status )
		solverErrmsg(status);

	mem_free(rowname[0]);
	mem_free(rowname);

	return status;
}//END addRow()

int addCol(LPptr lp, int nzcnt, double objx, int matbeg, intvec cmatind, vectorC cmatval, double bdu, double bdl, stringC *colname){
	int 	status, cmatbeg[1];
	double obj[1], lb[1], ub[1];

	cmatbeg[0] = matbeg;
	obj[0] = objx;
	ub[0] = bdu;
	lb[0] = bdl;

	status = CPXaddcols(env, lp, 1, nzcnt, obj, cmatbeg, cmatind, cmatval, lb, ub, NULL);
	if ( status )
		solverErrmsg(status);

	return status;
}//END addCol()

int removeRow(LPptr lp, int begin, int end){
	int status;

	status = CPXdelrows(env, lp, begin, end);
	if ( status )
		solverErrmsg(status);

	return status;
}//END removeRow

int freeProblem(LPptr lp) {
	int status;

	status = CPXfreeprob(env, &lp);
	if (status)
		solverErrmsg(status);

	return status;
}//END freeProb()
