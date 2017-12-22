/*
 * smps.h
 *
 *  Created on: Sep 29, 2015
 *      Author: Harsha Gangammanavar
 */

#ifndef SMPS_H_
#define SMPS_H_
#include "utils.h"

#undef INPUT_CHECK

typedef struct{
	int		type;			/* type of problem: LP, QP, MIP or MIQP */
	void	*lp;			/* problem pointer to be used by solver */
	string	name;			/* name of the problem */
	int		objsen;			/* sense of the objective: 1 for minimization and -1 for maximization */
	int		mac;			/* number of columns */
	int 	mar;			/* number of rows */
	int		numBin;			/* number of binary variables in the problem */
	int		numInt;			/* number of integer variables in the problem */
	int		numnz;			/* number of non-zero elements in constraint matrix */
	vector	objx;			/* objective function coefficients */
	vector	rhsx;			/* right-hand side */
	string	senx;			/* constraint sense */
	intvec	matbeg;			/* sparse matrix representation: column beginning */
	intvec	matcnt;			/* sparse matrix representation: number of non-zero entries in a column */
	intvec	matind;			/* sparse matrix representation: rows with non-zero entries */
	vector	matval;			/* sparse matrix representation: non-zero coefficients of the matrix */
	vector	bdl;			/* lower bound */
	vector	bdu;			/* upper bound */
	string	ctype;			/* type of decision variables: 'C' continuous, 'B' binary, 'I' general integer, 'S' semi-continuous, 'N' semi-integer */
	string	objname;		/* objective function name */
	int		rstorsz;		/* memory size for storing row names */
	string	*rname;			/* vector of row names */
	string	rstore;			/* row names string */
	int		cstorsz;		/* memory size for storing column names */
	string	*cname;			/* vector of column names */
	string	cstore;			/* column name string */
	int		macsz;			/* extended column size */
	int		marsz;			/* extended row size */
	int		matsz;			/* extended matrix size */
}oneProblem;

typedef struct {
	int		type;			      /* type of time file declaration, 0 for implicit and 1 for explicit */
	string	probName;		      /* name of the problem as read from time file */
	int		numStages;	       	  /* number of stages in the problem */
	string	*stgNames;		      /* unique strings to identify stages*/
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
	vector			meanEps;	  /* mean of noise process */
	vector			varEps;		  /* variance of noise process */
	vector			*eta;		  /* trend time series */
	vector			*sigma;		  /* seasonality time series */
	vector			*obs;		  /* historical time series (used for initialization) */
	vector			*eps;		  /* historical noise series (used for initilization) */
}armaType;

typedef struct {
	string	type;				/* type of stocType being used */
	BOOL	sim;				/* set to TRUE if an external simulator is used */
	int		numOmega; 			/* number of stochastic elements stored in structure */
	int 	numCipher; 			/* number of ints needed to encode an observation */
	int		numGroups;
	intvec	row; 				/* row number array in the original problem; -1 indicates objective function */
	intvec	col; 				/* column number array in the original problem; -1 indicates right-hand side */
	intvec	numVals;			/* number of realization for each random variable */
	vector	*vals; 				/* indexed array of discrete realizations of random variable */
	vector	*probs;				/* indexed array of probabilities associated with discrete realizations*/
	intvec	numPerGroup;
	intvec	groupBeg;
	vector	mean;         		/* mean of each rv */
	armaType *arma;
}stocType;

/* subroutines in smps.c */
int readFiles(string inputDir, string probName, oneProblem **orig, timeType **tim, stocType **stoc);
oneProblem *readCore(string inputDir, string probName);
timeType *readTime(string inputDir, string probName, oneProblem *orig);
stocType *readStoc(string inputDir, string probName, oneProblem *orig, timeType *tim);
int readIndep(FILE *fptr, string *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc);
int readBlocks(FILE *fptr, string *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc);
int readBlk(FILE *fptr, string *fields, oneProblem *orig, int maxOmegas, int maxVals, BOOL origRV, stocType *stoc);
int readARMA(FILE *fptr, string *fields, oneProblem *orig, stocType *stoc, int maxOmegas);
int readScenarios(FILE *fptr, string *fields, oneProblem *orig, timeType *tim, int maxOmegas, int maxVals, stocType *stoc);

void freeOneProblem(oneProblem *p);
void freeTimeType(timeType *tim);
void freeStocType(stocType *stoc);
void freeARMAtype(armaType *arma);

/* subroutines in rvgen.c */
int generateOmegaIdx(stocType *stoc, long long *seed);
void generateOmega(stocType *stoc, vector observ, long long *seed);
void generateBlocks(stocType *stoc, vector observ, int groupID, long long *seed);
void generateIndep(stocType *stoc, vector observ, int groupID, long long *seed);
int normal(vector mu, vector stdev, int numOmega, vector observ, long long *seed);
float scalit(float lower, float upper, long long *RUN_SEED);
float randUniform(long long *SEED);
int randInteger(long long *SEED, int iMax);
vector* setupSAA(stocType *stoc, long long *seed, int *numSamples);

#endif /* SMPS_H_ */
