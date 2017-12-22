/*
 * solver.h
 *
 *  Created on: Apr 20, 2014
 *      Author: gjharsha
 */
#ifndef MTSD_SOLVER_H_
#define MTSD_SOLVER_H_

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

int solveProblem(LPptr lp, string pname, int type, int *status);
int getProbType(LPptr lp);
double getObjective(LPptr lp, int type);
int getPrimal(LPptr lp, vector X, int length);
double getPrimalPoint(LPptr lp, int idx);
int getDual(LPptr lp, vector Pi, int length);
int getDualSlacks(LPptr lp, vector Dj, int length);
int getBasis(LPptr lp, intvec cstat, intvec rstat);
int getBinvC(LPptr lp, int col, vector a );
int changeCoef(LPptr lp, int row, int col, double val);
int changeObjx(LPptr lp, int cnt, intvec indices, vector values);
int changeRHS(LPptr lp, int cnt, intvec indices, vector values);
int changeBDS(LPptr lp, int cnt, intvec indices, string lu, vector bd);
int changeCol(LPptr lp, int column, vector coef, int start, int stop);
int changeCtype(LPptr lp, int cnt, intvec indices, string ctype);
int changeProbType(LPptr lp, int type);
int addRow(LPptr lp, int nzcnt, double inputRHS, char inputSense, int matbeg, intvec rmatind, vector rmatval);
int addCol(LPptr lp, int nzcnt, double objx, int cmatbeg, intvec cmatind, vector cmatval, double bdu, double bdl, string *colname);
int removeRow(LPptr lp, int begin, int end);

int createProblem(char *probname, LPptr *lp);
int readProblem(char *probpath, LPptr lp);
LPptr setupProblem(string name, int type, int numcols, int numrows, int objsense, vector objx, vector rhsx, string sense, intvec matbeg, intvec matcnt,
		intvec matind, vector matval, vector lb, vector ub, vector rngval, string *colname, string *rowname, string ctype);
int loadProblem(CPXLPptr lp, int numcols, int numrows, int objsense, vector objx, vector rhsx, string sense, intvec matbeg, intvec matcnt,
		intvec matind, vector matval, vector lb, vector ub, vector rngval);
int loadProbwNames(LPptr lp, int numcols, int numrows, int objsense, vector objx, vector rhsx, string sense, intvec matbeg, intvec matcnt,
		intvec matind, vector matval, vector lb, vector ub, vector rngval, string *colname, string *rowname);
int copyQPseparable(LPptr lp, double *qsepvec);
LPptr cloneProblem(LPptr origLp);
int writeProblem(LPptr lp, char *filename);

void openSolver();
void closeSolver();
int setIntParam(int paramname, int paramvalue);
void solverErrmsg(int status);
int changeSolverType();

int getProbName(LPptr lp, string probName, int len);
int getObjSen(LPptr lp);
int getNumRows(LPptr lp);
int getNumCols(LPptr lp);
int getNumBinary(LPptr lp);
int getCtype(LPptr lp, int start, int end, string ctype);
int getNumInt(LPptr lp);
int getNumnz(LPptr lp);
int getObjx(LPptr lp, int start, int end, vector obj);
int getRhsx(LPptr lp, int start, int end, vector rhs);
int getSense(LPptr lp, int start, int end, string sense);
int getCols(LPptr lp, int start, int end, int *cmatbeg, intvec cmatind, vector cmatval, int cmatspace);
int getLb(LPptr lp, int start, int end, vector lb);
int getUb(LPptr lp, int start, int end, vector ub);
int getObjName(LPptr lp, string objname);
int getCstoreSize(LPptr lp, int start, int end);
int getColName(LPptr lp, int start, int end, string *colname, string colnamestore, int csize);
int getRstoreSize(LPptr lp, int start, int end);
int getRowName(LPptr lp, int start, int end, string *rowname, string rownamestore, int rsize);
int getBinvC(LPptr lp, int col, vector a ) ;
int getBhead(LPptr lp, intvec head, vector x) ;
int binvArow(LPptr lp, int i, vector z) ;
int binvAcol(LPptr lp, int i, vector z);
int freeProblem(LPptr lp);

#endif /* MTSD_SOLVER_H_ */
