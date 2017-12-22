/*
 * prob.h
 *
 *  Created on: Sep 30, 2015
 *      Author: Harsha Gangammanavar
 */

#ifndef PROB_H_
#define PROB_H_

#undef DECOMPOSE_CHECK

/* Structure that holds various dimensions of the stage problem */
typedef struct {
	int		rows;			/* number of rows */
	int		cols;			/* number of columns */
	int		intCols;		/* number of integer variables */
	int		binCols;		/* number of binary variables */
	int		prevRows;		/* number of rows in previous stage, is set to 0 for root-stage (master) */
	int		prevCols;		/* number of columns in previous stage, is set to 0 for root-stage (master) */
	int     cntCcols;  		/* number of columns in T which have at least one non-zero element */
	int		numRV;			/* total number of random variables */
	int 	rvRowCnt;		/* number of rows effected by randomness */
	int 	rvColCnt;		/* number of columns effected by randomness */
    int 	rvaOmCnt;		/* number of RVs in dynamics noise vector a */
	int 	rvbOmCnt;		/* number of RVs in right-hand side */
    int 	rvcOmCnt;		/* number of RVs in state-cost coefficients */
    int 	rvdOmCnt;		/* number of RVs in stage-cost coefficients */
    int 	rvAOmCnt;		/* number of RVs in dynamics state-matrix A */
    int 	rvBOmCnt;		/* number of RVs in dynamics decision-matrix B */
	int 	rvCOmCnt;		/* number of RVs in transfer matrix C */
    int 	rvDOmCnt;		/* number of RVs in recourse matrix D */
}numType;

/* structure to hold coordinate information at each stage */
typedef struct {
	intvec	omegaRow;		/* list of all random variable rows */
	intvec	omegaCol;		/* list of all random variable columns */
	intvec	colsC;			/* list of columns in transfer matrix with at least non-zero element in them */
	intvec	rvCols;			/* list of all columns in transfer matrix with at least one random element */
	intvec	rvRows;			/* list of all rows with at least one random variable */
}coordType;

/* structure to hold exogenous information at each stage */
typedef struct {
	int		numRV;			/* total number of random variables */
	int		beg;			/* indicates the beginning of stage random variables in the entire observation vector */
	intvec	row;			/* row coordinate for random variable; -1 indicates cost coefficients */
	intvec	col;			/* column coordinates for random variable; -1 indicates right-hand side */
	vector	mean;			/* vector of mean values */
}omegastuff;

/* structure for the problem type:
 * c_t^\top x_t + \min  d_t^\top u_t + \expect{h_{t+}(s_{t+})}
 *                 s.t. D_t u_t = b_t - C_tx_t
 * where, x_{t+} = a_{t+} + A_{t+}x_t + B_{t+}u_t.
 */
typedef struct{
	string			name;			/* name of the problem */
	oneProblem		*sp;			/* structure with complete problem information */
	numType			*num;			/* structure which holds the problem dimensions */
	coordType		*coord;			/* structure which holds the necessary coordinates of the problem */
	omegastuff		*omegas;		/* exogenous information relevant to the problem */
	sparseVector	*aBar;			/* dynamics vector a_{t+} */
	sparseVector	*bBar;			/* right-hand side b_t */
	sparseVector	*cBar;			/* state cost coefficients c_t */
	sparseVector	*dBar;			/* objective function coefficients d_t */
	sparseMatrix	*Abar;			/* dynamics state matrix A_{t+} */
	sparseMatrix	*Bbar;			/* dynamics decision matrix B_{t+} */
	sparseMatrix	*Cbar;			/* transfer matrix C_t */
	sparseMatrix	*Dbar;			/* recourse matrix D_t */
	double			lb;				/* lower bounds on cost-to-go function */
}probType;

/* subroutines in prob.c */
probType **newProb(oneProblem *orig, stocType *stoc, timeType *tim, vector lb, double TOLERANCE);
vector meanProblem(oneProblem *orig, stocType *stoc);
vector calcLowerBound(oneProblem *orig, timeType *tim, stocType *stoc);
void freeProbType(probType **prob, int T);
void freeCoordType (coordType *coord);
void freeOmegastuff(omegastuff *omegas);
void printDecomposeSummary(timeType *tim, probType **prob);

#endif /* PROB_H_ */
