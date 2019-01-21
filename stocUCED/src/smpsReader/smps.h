/*
 * smps.h
 *
 *  Created on: Sep 29, 2015
 *      Author: Harsha Gangammanavar
 */

#ifndef SMPS_H_
#define SMPS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "../solverUtilities/utils.h"

#undef INPUT_CHECK

typedef struct{
	int		type;			/* type of problem: LP, QP, MIP or MIQP */
	void	*lp;			/* problem pointer to be used by solver */
	stringC	name;			/* name of the problem */
	int		objsen;			/* sense of the objective: 1 for minimization and -1 for maximization */
	int		mac;			/* number of columns */
	int 	mar;			/* number of rows */
	int		numBin;			/* number of binary variables in the problem */
	int		numInt;			/* number of integer variables in the problem */
	int		numnz;			/* number of non-zero elements in constraint matrix */
	vectorC	objx;			/* objective function coefficients */
	vectorC	rhsx;			/* right-hand side */
	stringC	senx;			/* constraint sense */
	intvec	matbeg;			/* sparse matrix representation: column beginning */
	intvec	matcnt;			/* sparse matrix representation: number of non-zero entries in a column */
	intvec	matind;			/* sparse matrix representation: rows with non-zero entries */
	vectorC	matval;			/* sparse matrix representation: non-zero coefficients of the matrix */
	vectorC	bdl;			/* lower bound */
	vectorC	bdu;			/* upper bound */
	stringC	ctype;			/* type of decision variables: 'C' continuous, 'B' binary, 'I' general integer, 'S' semi-continuous, 'N' semi-integer */
	stringC	objname;		/* objective function name */
	int		rstorsz;		/* memory size for storing row names */
	stringC	*rname;			/* vectorC of row names */
	stringC	rstore;			/* row names string */
	int		cstorsz;		/* memory size for storing column names */
	stringC	*cname;			/* vector of column names */
	stringC	cstore;			/* column name string */
	int		macsz;			/* extended column size */
	int		marsz;			/* extended row size */
	int		matsz;			/* extended matrix size */
}oneProblem;

typedef struct {
	int		type;			      /* type of time file declaration, 0 for implicit and 1 for explicit */
	stringC	probName;		      /* name of the problem as read from time file */
	int		numStages;	       	  /* number of stages in the problem */
	stringC	*stgNames;		      /* unique strings to identify stages*/
	intvec	row;			      /* a list of row names which mark the beginning of a new stage */
	intvec	col;			      /* a list of column names which mark the beginning of a new stage */
	int		numRows;		      /* used with explicit time file declaration only, set to numStages in implicit declaration */
	intvec	rowStg;			      /* used with explicit time file declaration only */
	int		numCols;		      /* used with explicit time file declaration only, set to numStages in implicit declaration */
	intvec	colStg;				  /* used with explicit time file declaration only */
}timeType;

typedef struct {
	int				p;			  /* autoregression order */
	int				q;			  /* moving-average order */
	int				N;			  /* dimension of time series */
	int				T;			  /* length of time series historical data */
	sparseMatrix	**AR;		  /* autoregression coefficients */
	sparseMatrix	**MA;		  /* moving-average coefficients */
	vectorC			meanEps;	  /* mean of noise process */
	vectorC			varEps;		  /* variance of noise process */
	vectorC			*eta;		  /* trend time series */
	vectorC			*sigma;		  /* seasonality time series */
	vectorC			*obs;		  /* historical time series (used for initialization) */
	vectorC			*eps;		  /* historical noise series (used for initilization) */
}armaType;

typedef struct {
	stringC	type;				/* type of stocType being used */
	BOOL	sim;				/* set to TRUE if an external simulator is used */
	int		numOmega; 			/* number of stochastic elements stored in structure */
	int 	numCipher; 			/* number of ints needed to encode an observation */
	int		numGroups;
	intvec	row; 				/* row number array in the original problem; -1 indicates objective function */
	intvec	col; 				/* column number array in the original problem; -1 indicates right-hand side */
	intvec	numVals;			/* number of realization for each random variable */
	vectorC	*vals; 				/* indexed array of discrete realizations of random variable */
	vectorC	*probs;				/* indexed array of probabilities associated with discrete realizations*/
	intvec	numPerGroup;
	intvec	groupBeg;
	vectorC	mean;         		/* mean of each rv */
	armaType *arma;
}stocType;

/* subroutines in smps.c */
int readFiles(stringC inputDir, stringC probName, oneProblem **orig, timeType **tim, stocType **stoc);
oneProblem *readCore(stringC inputDir, stringC probName);
timeType *readTime(stringC inputDir, stringC probName, oneProblem *orig);
stocType *readStoc(stringC inputDir, stringC probName, oneProblem *orig, timeType *tim);
int readIndep(FILE *fptr, stringC *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc);
int readBlocks(FILE *fptr, stringC *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc);
int readBlk(FILE *fptr, stringC *fields, oneProblem *orig, int maxOmegas, int maxVals, BOOL origRV, stocType *stoc);
int readARMA(FILE *fptr, stringC *fields, oneProblem *orig, stocType *stoc, int maxOmegas);
int readScenarios(FILE *fptr, stringC *fields, oneProblem *orig, timeType *tim, int maxOmegas, int maxVals, stocType *stoc);

void freeOneProblem(oneProblem *p);
void freeTimeType(timeType *tim);
void freeStocType(stocType *stoc);
void freeARMAtype(armaType *arma);

/* subroutines in rvgen.c */
int generateOmegaIdx(stocType *stoc, long long *seed);
void generateOmega(stocType *stoc, vectorC observ, long long *seed);
void generateBlocks(stocType *stoc, vectorC observ, int groupID, long long *seed);
void generateIndep(stocType *stoc, vectorC observ, int groupID, long long *seed);
int rndNormal(vectorC mu, vectorC stdev, int numOmega, vectorC observ, long long *seed);
float scalit(float lower, float upper, long long *RUN_SEED);
float randUniform(long long *SEED);
int randInteger(long long *SEED, int iMax);
vectorC* setupSAA(stocType *stoc, long long *seed, int *numSamples);

#ifdef __cplusplus
}
#endif

#endif /* SMPS_H_ */
