/*
 * solver.h
 *
 *  Created on: Apr 20, 2014
 *      Author: gjharsha
 */
#ifndef MTSD_SOLVER_H_
#define MTSD_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ilcplex/cplex.h"

#define		ENVptr			CPXENVptr
#define 	LPptr			CPXLPptr

#define		ON				CPX_ON
#define		OFF				CPX_OFF
#define 	INFBOUND    	CPX_INFBOUND

#define		PARAM_SCRIND	CPX_PARAM_SCRIND
#define		PARAM_SCAIND	CPX_PARAM_SCAIND
#define		PARAM_LPMETHOD	CPX_PARAM_LPMETHOD

#define		ALG_AUTOMATIC	CPX_ALG_AUTOMATIC
#define		ALG_PRIMAL		CPX_ALG_PRIMAL
#define		ALG_DUAL		CPX_ALG_DUAL
#define		ALG_NET			CPX_ALG_NET
#define		ALG_BARRIER		CPX_ALG_BARRIER
#define		ALG_SIFTING		CPX_ALG_SIFTING
#define		ALG_CONCURRENT	CPX_ALG_CONCURRENT

#define		STAT_OPTIMAL	CPX_STAT_OPTIMAL
#define		STAT_INFEASIBLE	CPX_STAT_INFEASIBLE

#define		AT_LOWER		CPX_AT_LOWER
#define		AT_UPPER		CPX_AT_UPPER

#define		MSGBUFSIZE		CPXMESSAGEBUFSIZE

#define		PROB_LP			CPXPROB_LP
#define		PROB_QP			CPXPROB_QP
#define		PROB_MILP		CPXPROB_MILP
#define		PROB_MIQP		CPXPROB_MIQP

#define		AT_LOWER        CPX_AT_LOWER
#define		BASIC           CPX_BASIC
#define		AT_UPPER        CPX_AT_UPPER
#define		FREE_SUPER      CPX_FREE_SUPER

#define		MIP_OPTIMAL		CPXMIP_OPTIMAL
#define		MIP_OPTIMAL_TOL	CPXMIP_OPTIMAL_TOL
#define		MIP_INFEASIBLE	CPXMIP_INFEASIBLE
#define     MIP_OPTIMAL_TOL CPXMIP_OPTIMAL_TOL

int solveProblem(LPptr lp, stringC pname, int type, int *status);
int getProbType(LPptr lp);
double getObjective(LPptr lp, int type);
int getPrimal(LPptr lp, vectorC X, int length);
double getPrimalPoint(LPptr lp, int idx);
int getDual(LPptr lp, vectorC Pi, int length);
int getDualSlacks(LPptr lp, vectorC Dj, int length);
int getBasis(LPptr lp, intvec cstat, intvec rstat);
int getBinvC(LPptr lp, int col, vectorC a );
int changeCoef(LPptr lp, int row, int col, double val);
int changeObjx(LPptr lp, int cnt, intvec indices, vectorC values);
int changeRHS(LPptr lp, int cnt, intvec indices, vectorC values);
int changeBDS(LPptr lp, int cnt, intvec indices, stringC lu, vectorC bd);
int changeCol(LPptr lp, int column, vectorC coef, int start, int stop);
int changeCtype(LPptr lp, int cnt, intvec indices, stringC ctype);
int changeProbType(LPptr lp, int type);
int addRow(LPptr lp, int nzcnt, double inputRHS, char inputSense, int matbeg, intvec rmatind, vectorC rmatval);
int addCol(LPptr lp, int nzcnt, double objx, int cmatbeg, intvec cmatind, vectorC cmatval, double bdu, double bdl, stringC *colname);
int removeRow(LPptr lp, int begin, int end);

int createProblem(char *probname, LPptr *lp);
int readProblem(char *probpath, LPptr lp);
LPptr setupProblem(stringC name, int type, int numcols, int numrows, int objsense, vectorC objx, vectorC rhsx, stringC sense, intvec matbeg, intvec matcnt,
		intvec matind, vectorC matval, vectorC lb, vectorC ub, vectorC rngval, stringC *colname, stringC *rowname, stringC ctype);
int loadProblem(CPXLPptr lp, int numcols, int numrows, int objsense, vectorC objx, vectorC rhsx, stringC sense, intvec matbeg, intvec matcnt,
		intvec matind, vectorC matval, vectorC lb, vectorC ub, vectorC rngval);
int loadProbwNames(LPptr lp, int numcols, int numrows, int objsense, vectorC objx, vectorC rhsx, stringC sense, intvec matbeg, intvec matcnt,
		intvec matind, vectorC matval, vectorC lb, vectorC ub, vectorC rngval, stringC *colname, stringC *rowname);
int copyQPseparable(LPptr lp, double *qsepvec);
LPptr cloneProblem(LPptr origLp);
int writeProblem(LPptr lp, char *filename);

void openSolver();
void closeSolver();
int setIntParam(int paramname, int paramvalue);
void solverErrmsg(int status);
int changeSolverType(int method);

int getProbName(LPptr lp, stringC probName, int len);
int getObjSen(LPptr lp);
int getNumRows(LPptr lp);
int getNumCols(LPptr lp);
int getNumBinary(LPptr lp);
int getCtype(LPptr lp, int start, int end, stringC ctype);
int getNumInt(LPptr lp);
int getNumnz(LPptr lp);
int getObjx(LPptr lp, int start, int end, vectorC obj);
int getRhsx(LPptr lp, int start, int end, vectorC rhs);
int getSense(LPptr lp, int start, int end, stringC sense);
int getCols(LPptr lp, int start, int end, int *cmatbeg, intvec cmatind, vectorC cmatval, int cmatspace);
int getLb(LPptr lp, int start, int end, vectorC lb);
int getUb(LPptr lp, int start, int end, vectorC ub);
int getObjName(LPptr lp, stringC objname);
int getCstoreSize(LPptr lp, int start, int end);
int getColName(LPptr lp, int start, int end, stringC *colname, stringC colnamestore, int csize);
int getRstoreSize(LPptr lp, int start, int end);
int getRowName(LPptr lp, int start, int end, stringC *rowname, stringC rownamestore, int rsize);
int getBinvC(LPptr lp, int col, vectorC a ) ;
int getBhead(LPptr lp, intvec head, vectorC x) ;
int binvArow(LPptr lp, int i, vectorC z) ;
int binvAcol(LPptr lp, int i, vectorC z);
int freeProblem(LPptr lp);

#ifdef __cplusplus
}
#endif

#endif /* MTSD_SOLVER_H_ */
