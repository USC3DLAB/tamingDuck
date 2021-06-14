/*
 * setup.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"

extern configType config;

int setupAlgo(oneProblem *orig, stocType *stoc, timeType *tim, probType ***prob, cellType **cell) {
	vectorC	meanSol, lb;
	int 	t;

	/* setup mean value problem which will act as reference for all future computations */
	meanSol = meanProblem(orig, stoc);
	if ( meanSol == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		return 1;
	}

	/* calculate lower bounds for each stage */
	lb = calcLowerBound(orig, tim, stoc);
	if ( lb == NULL )  {
		errMsg("setup", "setupAlgo", "failed to compute lower bounds on stage problem", 0);
		return 1;
	}

	/* decompose the problem into master and subproblem */
	(*prob) = newProb(orig, stoc, tim, lb, config.TOLERANCE);
	if ( (*prob) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		return 1;
	}

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(tim, (*prob));
#endif

	/* ensure that we have a linear programs at all stages */
	t = 0;
	while ( t < tim->numStages ) {
		if ( (*prob)[t++]->sp->type  != PROB_LP )
			printf("Warning :: Stage-%d problem is a mixed-integer program. Solving its linear relaxation.\n", t);
	}

	/* create the cells which will be used in the algorithms */
	(*cell) = newCell(stoc, (*prob), meanSol);
	if ( (*cell) == NULL ) {
		errMsg("setup", "setupAlgo", "failed to create the necessary cell structure", 0);
		return 1;
	}

	mem_free(meanSol); mem_free(lb);

	return 0;
}//END setupAlgo()

/* This function is used to create cells used in the algorithm */
cellType *newCell(stocType *stoc, probType **prob, vectorC xk) {
	cellType    *cell;
	int			length;

	/* allocate memory to all cells used in the algorithm. The first cell belongs to the master problem, while the rest correspond to each of the
	 * sub-agents in the problem.  */
	if (!(cell = (cellType *) mem_malloc(sizeof(cellType))) )
		errMsg("Memory allocation", "new_cell", "failed to allocate memory to cell",0);
	cell->master = cell->subprob = NULL;
	cell->candidX = cell->incumbX = NULL;
	cell->piS = cell->piM = cell->djM = NULL;
	cell->cuts = cell->fcuts = NULL;
	cell->lambda = NULL; cell->sigma = NULL; cell->delta = NULL; cell->omega = NULL;
	cell->pi_ratio = NULL;

	/* setup the master problem */
	cell->master = newMaster(prob[0], xk);
	if ( cell->master == NULL ) {
		errMsg("setup", "newCell", "failed to setup the master problem", 0);
		return NULL;
	}
	/* setup the subproblem */
	cell->subprob = newSubprob(prob[1]);

	/* -+-+-+-+-+-+-+-+-+-+-+ Allocating memory to other variables that belongs to master mcell +-+-+-+-+-+-+-+-+-+- */
	cell->k 	= 0;
	cell->LPcnt = 0;
	if (prob[0]->lb == 0)
		cell->lbType = TRIVIAL;
	else
		cell->lbType = NONTRIVIAL;
	cell->lb = prob[0]->lb;

	/* candidate solution and estimates */
	cell->candidX 			= duplicVector(xk, prob[0]->num->cols);
	cell->candidEst 		= prob[0]->lb + vXvSparse(cell->candidX, prob[0]->dBar);

	/* incumbent solution and estimates */
	if (config.MASTERTYPE == PROB_QP) {
		cell->incumbX   = duplicVector(xk, prob[0]->num->cols);
		cell->incumbEst = cell->candidEst;
		cell->quadScalar= config.MIN_QUAD_SCALAR;     						/* The quadratic scalar, 'sigma'*/
		cell->iCutIdx   = 0;
		cell->iCutUpdt  = 0;
		cell->incumbChg = CTRUE;
	}
	else {
		cell->incumbX   = NULL;
		cell->incumbEst = 0.0;
		cell->quadScalar= 0.0;
		cell->iCutIdx   = -1;
		cell->iCutUpdt  = -1;
		cell->incumbChg = CFALSE;
	}
	cell->gamma 			= 0.0;
	cell->normDk_1 			= 0.0;
	cell->normDk 			= 0.0;

	/* lower bounding approximations held in cuts structure */
	cell->maxCuts = config.CUT_MULT * prob[0]->num->cols + 3;
	cell->cuts 	  = newCuts(cell->maxCuts);
	cell->fcuts   = NULL;

	/* solution parts of the cell */
	if ( !(cell->piS = (vectorC) arr_alloc(prob[1]->num->rows + 1, double)) )
		errMsg("allocation", "newMaster", "cell->pi", 0);
	if ( !(cell->djM = (vectorC) arr_alloc(prob[0]->num->cols + 2, double)) )
		errMsg("allocation", "newMaster", "cell->di", 0);
	if ( !(cell->piM = (vectorC) arr_alloc(prob[0]->num->rows + cell->maxCuts + 1, double)) )
		errMsg("allocation", "newMaster", "cell->piM", 0);
	cell->mubBar = 0.0;

	/* stochastic elements */
	length = config.MAX_ITER + config.MAX_ITER / config.TAU + 1;
	cell->lambda = newLambda(length, 0, prob[1]->num->rvRowCnt);
	cell->sigma  = newSigma(length, prob[1]->num->rvColCnt, 0);
	cell->delta  = newDelta(length);
	cell->omega  = newOmega(config.MAX_ITER);

	cell->optFlag 			= CFALSE;
	cell->dualStableFlag 	= CFALSE;
	if ( !(cell->pi_ratio = (vectorC) arr_alloc(config.SCAN_LEN, double)) )
		errMsg("allocation", "newCell", "cell->pi_ratio", 0);

	cell->feasCnt 			= 0;
	cell->infeasIncumb 		= CFALSE;

	/* construct the QP using the current incumbent */
	if ( config.MASTERTYPE == PROB_QP && cell->incumbChg == CTRUE) {
		/* update the right-hand side and the bounds with incumbent solution */
		if ( changeQPrhs(prob[0], cell) ) {
			errMsg("setup", "newCell", "failed to change the right-hand side after incumbent change", 0);
			return NULL;
		}
		if ( changeQPbds(CPXLPptr(cell->master->lp), prob[0]->num->cols, prob[0]->sp->bdl, prob[0]->sp->bdu, cell->incumbX) ) {
			errMsg("setup", "newCell", "failed to change the bounds after incumbent update", 0);
			return NULL;
		}
		cell->incumbChg = CFALSE;

		/* change the proximal term */
		if ( changeQPproximal(CPXLPptr(cell->master->lp), prob[0]->num->cols, config.MIN_QUAD_SCALAR) ) {
			errMsg("setup", "newCell", "failed to add the proximal term to QP", 0);
			return NULL;
		}
	}

	return cell;
}//END newCell()

void freeCellType(cellType *cell) {

	if ( cell ) {
		if (cell->master) freeOneProblem(cell->master);
		if (cell->candidX) mem_free(cell->candidX);
		if (cell->incumbX) mem_free(cell->incumbX);
		if (cell->piS) mem_free(cell->piS);
		if (cell->piM) mem_free(cell->piM);
		if (cell->djM) mem_free(cell->djM);
		if (cell->cuts) freeCutsType(cell->cuts);
		if (cell->fcuts) freeCutsType(cell->fcuts);
		if (cell->delta) freeDeltaType(cell->delta, cell->lambda->cnt, cell->omega->cnt);
		if (cell->omega) freeOmegaType(cell->omega);
		if (cell->lambda) freeLambdaType(cell->lambda);
		if (cell->sigma) freeSigmaType(cell->sigma);
		if (cell->pi_ratio) mem_free(cell->pi_ratio);
		mem_free(cell);
	}

}//END freeCellType()
