/*
 * twoSD.h
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#ifndef TWOSD_H_
#define TWOSD_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

#define TRIVIAL 0
#define NONTRIVIAL 1
#define INF	DBL_MAX

#undef STOCH_CHECK
#undef ALGO_CHECK

typedef struct{
	long long RUN_SEED;			/* seed used during optimization */
	double 	TOLERANCE; 			/* for zero identity test */
	int		MIN_ITER;			/* minimum number of iterations */
	int		MAX_ITER;			/* maximum number of iterations */
	int		MASTERTYPE;			/* type of master problem */
	int		CUT_MULT;			/* Determines the number of cuts to be used for approximate */
	int		TAU;				/* Frequency at which the incumbent is updated */
	double	MIN_QUAD_SCALAR;	/* Minimum value for regularizing parameter */
	double 	MAX_QUAD_SCALAR;	/* Maximum value for regularizing parameter */
	double	R1;
	double	R2;
	double	R3;
	int		PI_EVAL_START;
	int		PI_CYCLE;
	int		SCAN_LEN;
	int		EVAL_FLAG;
	long long EVAL_SEED;
	int		EVAL_MIN_ITER;
	double	EVAL_ERROR;
	double  PRE_EPSILON;
	double	EPSILON;
	int		BOOTSTRAP_REP;		/* Number of boot-strap replications in full optimality test */
	double	PERCENT_PASS;		/* percentage of bootstrap replications need to be satisfied */
}configType;

typedef struct {
	double  alpha;                  /* scalar value for the righ-hand side */
	vectorC  beta;                   /* coefficients of the master problems's primal variables */
	int 	cutObs;					/* number of samples on which the given cut was based */
	int 	omegaCnt;				/* number of *distinct* observations on which the cut is based (this is also the length of istar) */
	intvec	iStar;					/* indices of maximal pi for each distint observation */
	BOOL	isIncumb;				/* indicates if the cut is an incumbent cut */
	double 	alphaIncumb;			/* right-hand side when using QP master, this is useful for quick updates */
	int 	slackCnt;				/* number of times a cut has been slack, used in deciding when the cut needs to be dropped */
	int 	rowNum;					/* row number for master problem in solver */
}oneCut;

typedef struct {
	int     cnt;                    /* number of cuts */
	oneCut  **vals;
}cutsType;

/* To save time and space, Pi x b and Pi x C are calculated as soon as possible and stored in structures like sigma and delta.  Toward
 * this end, pixbCType represents a single calculation of pi X b (which is a scalar) and pi X C (which is a vectorC).*/
typedef struct{
	double 	pib;
	vectorC 	piC;
} pixbCType;

/* The lambda structure stores some of the dual variable values from every distinct dual vectorC obtained during the program.  Each vectorC contains
 * only those dual variables whose corresponding rows in the subproblem constraint matrix contain random elements.  _val_ is an array of
 * these dual vectors (thus it is 2-D). _row_ gives the corresponding row number for a given dual variable in _val_.  _cnt_ represents
 * the number of dual vectors currently stored in lambda. */
typedef struct {
	int 	cnt;
	vectorC 	*vals;
} lambdaType;

/* The sigma matrix contains the values of Pi x bBar and Pi x Cbar  for all values of pi obtained so far (note it does not depend on
 * observations of omega).  _col_ gives the column number of each non-zero element in pi X Cbar.  _val_ is an array of values
 * for pi X bBar and pi X Cbar, one entry for each pi.  Note that values  which are always zero (because Rbar or Cbar is zero there) are not
 * stored.  The _lamb_ array is the same size as the _val_ array, and for each element in _val_ the corresponding element in _lamb_ references
 * the dual vectorC in lambda that was used to calculate that entry in sigma. */
typedef struct {
	int 		cnt;
	pixbCType 	*vals;
	intvec		lambdaIdx;
	intvec		ck; 				/* record the iteration # when sigma was created */
} sigmaType;

/* When calculating istar for a cut, it is useful to have two separate references into the sigma and delta structures, since each dual vector
 * is stored in two places -- part in sigma and part in delta.  The final entry in cut->istar[] will just be the _sigma_ field of this structure. */
typedef struct {
	int delta;
	int sigma;
} iType;

/* The delta matrix contains the values of lambda_pi X bOmega and lambda_pi X Comega for all values of pi and all observations of omega.
 * _col_ gives the column number of the each non-zero element in the multiplication of lambda_pi X Comega (the same elements are non-zero
 * each time).  _val_ is an array of vectors of (lambda_pi X bOmega, lambda_pi X Comega) pairs of calculations.  A row in _val_ corresponds
 * to a distinct dual vector, and a column in _val_ corresponds to a distinct observation of omega.  Thus, every pi-omega combination is
 * represented here, and the size of the delta matrix can be determine from lambda->cnt and omega->cnt.
 * Note that when elements of omega get dropped, vacant columns appear in delta.  This is ok, but be sure to loop carefully! */
typedef struct {
	pixbCType 	**vals;
} deltaType;

/**************************************************************************\
 ** Omega stores the set of observations which have been made so far.
 **
 **   Each observation consists of a vectorC of realizations of random
 ** variables with discrete distributions.  Since every distribution is
 ** discrete, an observation is just stored as a vectorC of indices into a
 ** distribution array containing the possible values.  _idx_ is an array
 ** of such vectors.  Each rv occurs in the R vectorC or T matrix which
 ** (along with the candidate X) make up the rhs of the subproblem.
 **
 **   The _row_ and _col_ arrays give the coordinates in R and T of each rv
 ** realization in a vector.  If _col_ is zero for a given entry, then the
 ** realization comes from R; otherwise, it comes from T.  The field _RT_
 ** represents a vectorC of actual realizations (as opposed to indices) of
 ** omega for both the Romega and Tomega structures (only one observation's
 ** worth of omega).
 **
 **   The _weight_ field specifies the number of times a particular outcome
 ** has been observed (it starts at 1 when the outcome is first generated,
 ** and increments every time the same outcome is observed again).  _cnt_
 ** just specifies the number of distinct outcomes which have been observed
 ** and stored in the omega structure.
 **
 **   If the problem gets too large, some unexciting omegas may be dropped.
 ** So, _filter_ (same length as _weight_) describes which vectors in the
 ** _idx_ array are actually filled with observations.  _next_ references
 ** the place to start in the _filter_ array when trying to find the next
 ** available position in _idx_.  _most_ represents the number of
 ** elements needed to store all elements in the _idx_ array, from
 ** the first to the last.  (It is the greatest index at which an
 ** observation is stored in the _idx_ array, plus one).   Finally,
 ** _last_ indicates the index of omega from which we started dropping
 ** omegas last time -- so we only "need" to go from omega.most down
 ** to omega.last when dropping them this time (we'll miss some).
 \**************************************************************************/
typedef struct {
	int 	cnt;
	intvec	weight;                 /* number of times that an omega is observed */
	vectorC	*vals;
} omegaType;

typedef struct{
    double		iterSolTime;        /* time of solving prob in each iter */
    double		totSolTime;         /* time of solving prob(cumulated) */
    double		iterCutGenTime;         /* time of generating cuts for each time */
    double  	totCutGenTime;      /* time of generating cuts(cumulated) */
    double		iterTime;        /* time for each iter (solveAgent) */
    double      totTime;     /* total time of doing a iteration for subprob */
//    vectorC  	masterIterTime;     /* time for each iter (within the while loop) */
//    double      totMasterIterTime;  /* time of cumulated iteration time */
}runTimeType;

typedef struct {
	int         k;                  /* number of iterations */
	int 		LPcnt; 				/* the number of LPs solved. */
    double		lb;					/* lower bound on cell objective function */
    int			lbType;				/* type of lower bound being used TRIVIAL if 0, else NONTRIVIAL */

    oneProblem  *master;            /* store master information */
	oneProblem 	*subprob;			/* store subproblem information */

	vectorC      candidX;            /* primal solution of the master problem */
	double      candidEst;          /* objective value master problem */

	vectorC      incumbX;			/* incumbent master solution */
	double      incumbEst;			/* estimate at incumbent solution */
	double 		quadScalar; 		/* the proximal parameter/quadratic scalar 'sigma' */
	BOOL        incumbChg;			/* set to be true if the incumbent solution has changed in an iteration */
	int         iCutIdx;			/* index of incumbent cut in cell->cuts structure */
	int         iCutUpdt;			/* iteration number when incumbent cut is updated */
	double      gamma;				/* improvement in objective function value */
	double      normDk_1;			/* (\Delta x^{k-1})^2 */
	double      normDk;				/* (\Delta x^k)^2 */

	vectorC      piS;                 /* subproblem dual information */
	double      mubBar;				/* dual slack information for subproblem */
	vectorC 		piM;				/* master dual information */
	vectorC      djM;                /* master reduced cost vectorC */

    int      	maxCuts;            /* maximum number of cuts to be used*/
	cutsType    *cuts;              /* optimality cuts */
	cutsType    *fcuts;             /* feasibility cuts */

	lambdaType 	*lambda;			/* holds dual solutions corresponding to rows effected by randomness */
	sigmaType 	*sigma;				/* holds $\pi \times \bar{b}$ and $\pi \times \bar{C} $ values */
	deltaType   *delta;				/* calculations based on realization and dual solutions observed */
	omegaType 	*omega;				/* all realizations observed during the algorithm */

    BOOL        optFlag;
	vectorC      pi_ratio;
    BOOL        dualStableFlag; 	/* indicates if dual variables are stable */

	int			feasCnt;			/* keeps track of the number of times infeasible candidate solution was encountered */
	BOOL		infeasIncumb;		/* indicates if the incumbent solution is infeasbible */
}cellType;

/* twoSD.c */
int readConfig();

/* algo.c */
int algo(oneProblem *orig, timeType *tim, stocType *stoc, stringC probName);
int solveCell(stocType *stoc, probType **prob, cellType *cell, stringC probName);
void writeStatistic(FILE **soln, probType *prob, cellType *cell, stringC probName);
void cleanupAlgo(probType **prob, cellType *cell, int T);

/* setup.c */
int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell);
cellType *newCell(stocType *stoc, probType **prob, vectorC xk);
void freeCellType(cellType *cell);

/* master.c */
int solveQPMaster(numType *num, sparseVector *dBar, cellType *cell, int IniRow, double lb);
int addCut2Master(cellType *cell, oneCut *cut, int lenX, double lb);
int constructQP(probType *prob, LPptr lp, vectorC incumbX);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb);
int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb);
int changeQPproximal(LPptr lp, int numCols, double sigma);
int changeQPrhs(probType *prob, cellType *cell);
int changeQPbds(LPptr lp, int numCols, vectorC bdl, vectorC bdu, vectorC xk);
oneProblem *newMaster(probType *prob, vectorC xk);

/* cuts.c */
int formSDCut(probType *prob, cellType *cell, vectorC Xvect, int omegaIdx, BOOL newOmegaFlag);
oneCut *SDCut(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, omegaType *omega, vectorC Xvect, int numSamples,
		BOOL *dualStableFlag, vectorC pi_ratio, double lb);
iType computeIstar(numType *num, coordType *coord, sigmaType *sigma, deltaType *delta, vectorC Xvect, vectorC PiCbarX, int obs,
		int ictr, BOOL pi_eval, double *argmax);
iType compute_new_istar(int obs, oneCut *cut, sigmaType *sigma, deltaType *delta, vectorC Xvect, numType *num, coordType *coord,
		vectorC PiCbarX, double *argmax, int ictr);
oneCut *newCut(int numX, int numIstar, int numSamples);
cutsType *newCuts(int maxCuts);
int reduceCuts(cellType *cell, vectorC candidX, vectorC pi, int betaLen, double lb);
int dropCut(cellType *cell, int cutIdx);
double calc_var(double *x, double *mean_value, double *stdev_value, int batch_size);
void print_cut(cutsType *cuts, numType *num, int idx);
void freeOneCut(oneCut *cut);
void freeCutsType(cutsType *cuts);
double calc_var(double *x, double *mean_value, double *stdev_value, int batch_size);

/* subprob.c */
int solveSubprob(probType *prob, cellType *cell, vectorC Xvect, int omegaIdx, BOOL newOmegaFlag);
vectorC computeRHS(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, vectorC X, vectorC obs);
void chgRHSwSoln(sparseVector *bBar, sparseMatrix *Cbar, vectorC rhs, vectorC X) ;
int chgRHSwObserv(LPptr lp, numType *num, coordType *coord, vectorC observ, vectorC spRHS, vectorC X);
oneProblem *newSubprob(probType *subprob);

/* stocUpdate.c */
int stochasticUpdates(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *Cbar, lambdaType *lambda, sigmaType *sigma,
                       deltaType *delta, omegaType *omega, BOOL newOmegaFlag, int omegaIdx, int maxIter, int iter, vectorC pi, double mubBar);
void calcDeltaCol(numType *num, coordType *coord, lambdaType *lambda, vectorC observ, int omegaIdx, deltaType *delta);
int calcLambda(numType *num, coordType *coord, vectorC Pi, lambdaType *lambda, BOOL *newLambdaFlag);
int calcSigma(numType *num, coordType *coord, sparseVector *bBar, sparseMatrix *CBar, vectorC pi, double mubBar,
              int idxLambda, BOOL newLambdaFlag, int iter, sigmaType *sigma, BOOL *newSigmaFlag);
int calcDeltaRow(int maxIter, numType *num, coordType *coord, omegaType *omega, lambdaType *lambda, int lambdaIdx, deltaType *delta);
int calcOmega(vectorC observ, int begin, int end, omegaType *omega, BOOL *newOmegaFlag);
int computeMU(LPptr lp, int numCols, double *mubBar);
lambdaType *newLambda(int num_iter, int numLambda, int numRVrows);
sigmaType *newSigma(int numIter, int numNzCols, int numPi);
deltaType *newDelta(int numIter);
omegaType *newOmega(int numIter);
void freeLambdaType(lambdaType *lambda);
void freeSigmaType(sigmaType *sigma);
void freeOmegaType(omegaType *omega);
void freeDeltaType (deltaType *delta, int lambdaCnt, int omegaCnt);

/* soln.c */
int checkImprovement(probType *prob, cellType *cell, int candidCut);
int replaceIncumbent(probType *prob, cellType *cell, double candidEst);
double maxCutHeight(cutsType *cuts, int currIter, vectorC xk, int betaLen, double lb);
double cutHeight(oneCut *cut, int currIter, vectorC xk, int betaLen, double lb);

/* optimal.c */
BOOL optimal(probType **prob, cellType *cell);
BOOL preTest(cellType *cell);
BOOL fullTest(probType **prob, cellType *cell);
cutsType *chooseCuts(cutsType *cuts, vectorC pi, int lenX);
void reformCuts(sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord, cutsType *gCuts, int *observ, int k, int lbType, int lb, int lenX);
double calcBootstrpLB(probType *prob, vectorC incumbX, vectorC piM, vectorC djM, int currIter, double quadScalar, cutsType *cuts);
void empiricalDistribution(omegaType *omega, int *cdf);
void resampleOmega(intvec cdf, intvec observ, int numSamples);

/* evaluate.c */
int evaluate(FILE **soln, stocType *stoc, probType **prob, cellType *cell, vectorC Xvect);

#ifdef __cplusplus
}
#endif

#endif /* TWOSD_H_ */
