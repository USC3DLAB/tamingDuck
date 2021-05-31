/*
 * soln.c
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
extern FILE* file;

/***********************************************************************\
 ** This function determines whether the "stagewise descent property" is
 ** satisified.  If the current approximation of f_k gives a lower difference
 ** between the candidate and incumbent x than the previous approximation
 ** gave, then the incumbent x is updated to the candidate x, and the
 ** reference to the incumbent cut is updated as well.  The function returns
 ** TRUE if the incumbent was updated; FALSE otherwise.
 \***********************************************************************/
int checkImprovement(probType *prob, cellType *cell, int candidCut) {
	double  candidEst;

	/* Calculate height at new candidate x with newest cut included */
	candidEst = vXvSparse(cell->candidX, prob->dBar) + maxCutHeight(cell->cuts, cell->k, cell->candidX, prob->num->cols, cell->lb);;
	cell->incumbEst = vXvSparse(cell->incumbX, prob->dBar) + maxCutHeight(cell->cuts, cell->k, cell->incumbX, prob->num->cols, cell->lb);

#ifdef SOL_CHECK
	printf("AggcandidEst =%lf, AggIncumEst =%lf\n",AggcandidEst, cell->incumbEst);
#endif

	/* If we see considerable improvement, then change the incumbent */
	if ((candidEst - cell->incumbEst) < (config.R1 * cell->gamma)) {
		/* when we find an improvement, then we need to replace the incumbent x with candidate x */
		if ( replaceIncumbent(prob, cell, candidEst) ) {
			errMsg("algorithm", "checkImprovement", "failed to replace incumbent solution with candidate", 0);
			return 1;
		}
		cell->iCutIdx = candidCut;
		cell->incumbChg = CFALSE;
		fprintf(file, "+"); fflush(file);
	}
	else {
		/* Update quad_scalar when no incumbent is found. */
		cell->quadScalar = min(config.MAX_QUAD_SCALAR, cell->quadScalar / config.R2);
		cell->normDk_1 = cell->normDk;
	}

	if ( changeQPproximal(CPXLPptr(cell->master->lp), prob->num->cols, cell->quadScalar) ) {
		errMsg("setup", "newCell", "failed to add the proximal term to QP", 0);
		return 1;
	}

	return 0;
}//END checkImprovement()

int replaceIncumbent(probType *prob, cellType *cell, double candidEst) {

	/* replace the incumbent solution with the candidate solution */
	copyVector(cell->candidX, cell->incumbX, prob->num->cols, CTRUE);
	cell->incumbEst = candidEst;

	/* update the proximal parameter based on estimated improvement */
	if ( cell->normDk > config.TOLERANCE )
		if ( cell->normDk >= config.R3 * cell->normDk_1 ) {
			cell->quadScalar *= config.R2 * config.R3 * cell->normDk_1/ cell->normDk;
			cell->quadScalar  = min(config.MAX_QUAD_SCALAR, cell->quadScalar);
			cell->quadScalar = max(config.MIN_QUAD_SCALAR, cell->quadScalar);
		}

	/* update the right-hand side and the bounds with new incumbent solution */
	if ( changeQPrhs(prob, cell) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the right-hand side after incumbent change", 0);
		return 1;
	}
	if ( changeQPbds(CPXLPptr(cell->master->lp), prob->num->cols, prob->sp->bdl, prob->sp->bdu, cell->incumbX) ) {
		errMsg("algorithm", "replaceIncumbent", "failed to change the bounds after incumbent update", 0);
		return 1;
	}

	/* update the candidate cut as the new incumbent cut */
	cell->iCutUpdt = cell->k;
	cell->incumbChg = CTRUE;

	/* keep the two norm of solution*/
	cell->normDk_1 = cell->normDk;
	/* Since incumbent solution is now replaced by a candidate, we assume it is feasible now */
	cell->infeasIncumb = CFALSE;
	/* gamma needs to be reset to 0 since there's no difference between candidate and incumbent*/
	cell->gamma = 0.0;

	return 0;
}//END replaceIncumbent()

/* This function loops through a set of cuts and find the highest cut height at the specified position x */
double maxCutHeight(cutsType *cuts, int currIter, vectorC xk, int betaLen, double lb) {
	double Sm = -INF, ht = 0.0;
	int cnt;

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		ht = cutHeight(cuts->vals[cnt], currIter, xk, betaLen, lb);
		if (Sm < ht) {
			Sm = ht;
		}
	}

	return Sm;
}//END maxCutHeight

/* This function calculates and returns the height of a given cut at a given X.  It includes the k/(k-1) update, but does not include
 * the coefficients due to the cell. */
double cutHeight(oneCut *cut, int currIter, vectorC xk, int betaLen, double lb) {
	double height;
	double t_over_k = ((double) cut->cutObs / (double) currIter);

	/* A cut is calculated as alpha - beta x X */
	height = cut->alpha - vXv(cut->beta, xk, NULL, betaLen);

	/* Weight cut based on number of observations used to form it */
	height *= t_over_k;

	/* Updated for optimality cut height*/
	height += (1 - t_over_k) * lb;

	return height;
}//END cutHeight()

