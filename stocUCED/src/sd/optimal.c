/*
 * optimal.c
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
extern FILE *file;

/* This function determines whether or not the current incumbent solution is considered to be optimal. Optimality is guarenteed if the
 * following criteria are satisfied:
 * 		0. Minimum number of iterations have been completed.
 * 		1. Dual solution set has stabilized.
 * 		2. If dual solution set is stable, the pre-test checks for "convergence" of objective function estimate.
 * 		3. Full test is based on boot-strapping, and checks the gap between primal (upper) and dual (lower) values.
 * The pre-test is performed only after the dual solution set has stabilized, and the full test is performed only if the pre-test is successful. */
BOOL optimal(probType **prob, cellType *cell) {

	/* ensure that the minimum number of iterations have been completed */
	if (cell->k > config.MIN_ITER && cell->dualStableFlag ) {
		/* perform the pre-test */
		if ( preTest(cell) ) {
			if (fullTest(prob, cell)) {
				/* full test satisfied */
				cell->optFlag = TRUE;
				fprintf (file, ">"); fflush(file);
				return TRUE;
			}
			else {
				fprintf(file, ">"); fflush(file);
			}

		}
	}

	return FALSE;
}//optimal()


/* Because checking optimality is an arduous task, we first do a pre-check to determine if the full test is worthwhile. This function
 * determines whether the height at the candidate is close enough to the height at the incumbent to warrant an optimality test. */
BOOL preTest(cellType *cell) {

	/* The candidate must be within some small percentage of incumbent cut */
	/* rare situation for cell->candid_est < 0 and cell->incumb_est > 0 */
	/* Note: cell->candidEst and cell->incumbEst could be 0 */
	if (cell->candidEst >= 0){
		cell->optFlag = (cell->candidEst >= (1 - config.PRE_EPSILON) * cell->incumbEst);
	}
	else
		cell->optFlag = (cell->candidEst > (1 + config.PRE_EPSILON) * cell->incumbEst);

	return cell->optFlag;

}//END preTest()

/* This function performs a complete statistical test of optimality. First, it selects cuts whose height at the incumbent is "close" to
 * incumbent cut's height.  Then, it performs M resamplings of the observations in omega, and reforms the selected cuts with respect
 * to these observations (as if *they* were observed instead of the actual omega).  For each of the M resamplings, a master program
 * containing the reformed cuts is solved, and if almost all of the solutions to these master programs.
 * Multi-cut version selects "good" cuts within each agent, and reforms theses cuts. Then the cut value are aggregated.
 * Single-cut version first selects "good" cuts in the master prob. Then by using the information of these "good" cuts to find the
 * corresponding cuts in each agent. After reforming these agent cuts, aggregated these reformed cuts as one single cut and add it to the
 * master problem */
BOOL fullTest(probType **prob, cellType *cell) {
	cutsType *gCuts;
	intvec  cdf, observ;
	double  est, ht, LB=0.0;
	int 	numPass = 0, rep, j;

	/* (a) choose good cuts */
	gCuts = chooseCuts(cell->cuts, cell->piM, prob[0]->num->cols);

	/* (b) calculate empirical distribution of omegas */
	if ( !(cdf = (intvec) arr_alloc(cell->omega->cnt+1, int)) )
		errMsg("allocation", "fullTest", "failed to allocate memory to cdf",0);
	if ( !(observ = (intvec) arr_alloc(cell->k, int)))
		errMsg("allocation", "fullTest", "resampled observations", 0);

	empiricalDistribution(cell->omega, cdf);

	for (rep = 0; rep < config.BOOTSTRAP_REP; rep++) {
		/* (c) resample from the set of observations */
		resampleOmega(cdf, observ, cell->k-1);

		/* (d) reform the good cuts by plugging in the omegas */
		reformCuts(cell->sigma, cell->delta, cell->omega, prob[1]->num, prob[1]->coord, gCuts, observ, cell->k-1, cell->lbType, prob[0]->lb, prob[0]->num->cols);

		/* (e) find out the best reformed cut estimate at the incumbent solution */
		est = gCuts->vals[0]->alpha - vXv(gCuts->vals[0]->beta, cell->incumbX, NULL, prob[0]->num->cols);
		for (j = 1; j < gCuts->cnt; j++) {
			ht = gCuts->vals[j]->alpha - vXv(gCuts->vals[j]->beta, cell->incumbX, NULL, prob[0]->num->cols);
			if ( est < ht)
				est = ht;
		}

		/* (f) Solve the master with reformed "good cuts" (all previous cuts are dropped) to obtain a lowe bound. In QP approach, we don't include the incumb_x * c in estimate */
		if (config.MASTERTYPE == 1)
			est += vXvSparse(cell->incumbX, prob[0]->dBar);
		else
			LB = calcBootstrpLB(prob[0], cell->incumbX, cell->piM, cell->djM, cell->k, cell->quadScalar, gCuts);

#if 0
		printf("\niter = %d, replication = %d, UB = %f, LB = %f, Gap = %lf", cell->k, rep, est, LB, DBL_ABS((est - LB) / cell->incumbEst));
#endif

		/* (g) compare the normalized difference between estimate and the lower bound. If the problem is a QP problem, we don't need add the constant term c^T x \hat{x} */
		if (DBL_ABS((est - LB) / cell->incumbEst) <= config.EPSILON)
			numPass++;

		/* (h) check No. of fails. skip out of the loop if there's no hope of meeting the condition */
		if ( rep + 1 - numPass >= (1 - config.PERCENT_PASS) * config.BOOTSTRAP_REP) {
			/* The bootstrap test has failed */
			mem_free(cdf); mem_free(observ); freeCutsType(gCuts);
			return FALSE;
		}
	}//END replication loop

	mem_free(cdf); mem_free(observ);
	freeCutsType(gCuts);
	return TRUE;

}//END full_test()

/* This function selects all cuts whose height at the incumbent solution is close to the height of the incumbent cut. These cuts together
 * are likely to provide good approximations of the recourse function at incumbent solution, when they are reformed with new observations.
 * The function returns a new cut structure which contains room for cuts to be reformed. Only the _istar_ and _cut_obs_ fields of
 * each cut have been initialized. */
cutsType *chooseCuts(cutsType *cuts, vectorC pi, int lenX) {
	cutsType *gCuts;
	int cnt;

	gCuts = newCuts(cuts->cnt);

	for ( cnt = 0; cnt < cuts->cnt; cnt++ ) {
		if (pi[cuts->vals[cnt]->rowNum + 1] > config.TOLERANCE) {
			gCuts->vals[gCuts->cnt] = newCut(lenX, cuts->vals[cnt]->omegaCnt, cuts->vals[cnt]->cutObs);
			copyIntvec(cuts->vals[cnt]->iStar, gCuts->vals[gCuts->cnt]->iStar, cuts->vals[cnt]->omegaCnt);
			gCuts->vals[gCuts->cnt]->rowNum = cuts->vals[cnt]->rowNum;
			gCuts->cnt++;
		}
	}

	return gCuts;
}//END choose_cuts

/* This function forms an empirical distribution on the observations stored in omega, and calculates an integer cdf to represent the distribution.
 * An observation which has been seen n times will have n times the probability of being chosen as an observation seen only once. */
void empiricalDistribution(omegaType *omega, intvec cdf) {
	int cnt;

	/* Calculate an integer cdf distribution for observations */
	cdf[0] = omega->weight[0];
	for (cnt = 1; cnt < omega->cnt; cnt++)
		cdf[cnt] = cdf[cnt - 1] + omega->weight[cnt];

}//END empirical_distrib

/* This function randomly selects a new set of observations from the old set of observations stored in omega.  Entries in omega which have been observed
 * multiple times have a proportionally higher chance of being selected for the new set.  The function fills an array, assumed to be of a size equal to the
 * number of iterations, with the new set of observations. */
void resampleOmega(intvec cdf, intvec observ, int numSamples) {
	int cnt, obs;
	int sample;

	/* Choose k observations according to cdf (k = number of iterations) */
	for (obs = 0; obs < numSamples; obs++) {
		sample = randInteger(&config.EVAL_SEED, numSamples);
		for (cnt = 0; sample > cdf[cnt]; cnt++)
			/* Loop until sample falls below cdf */;
		observ[obs] = cnt;
	}
}//END resampleOmega

/* This function will calculate a new set of cuts based on the observations of omega passed in as _observ_, and the istar's which have already been stored in
 * the _istar_ field of each cut. If an istar field does not exist for a given observation, then a value of zero is averaged into the calculation of alpha & beta. */
void reformCuts(sigmaType *sigma, deltaType *delta, omegaType *omega, numType *num, coordType *coord, cutsType *gCuts, int *observ, int k, int lbType, int lb, int lenX) {
	int cnt, obs, idx, count;
	iType iStar;

	/* Loop through all the cuts and reform them */
	for (cnt = 0; cnt < gCuts->cnt; cnt++) {
		/* Begin with cut coefficients of zero */
		for (idx = 0; idx <= lenX; idx++)
			gCuts->vals[cnt]->beta[idx] = 0.0;
		gCuts->vals[cnt]->alpha = 0.0;

		count = 0;
		/* Reform this cut based on resampled observations */
		for (obs = 0; obs < k; obs++) {
			/* Only sum values if the cut has an istar for this observation */
			if (observ[obs] < gCuts->vals[cnt]->omegaCnt) {
				iStar.sigma = gCuts->vals[cnt]->iStar[observ[obs]];
				iStar.delta = sigma->lambdaIdx[iStar.sigma];

				gCuts->vals[cnt]->alpha += sigma->vals[iStar.sigma].pib + delta->vals[iStar.delta][observ[obs]].pib;

				for (idx = 1; idx <= num->cntCcols; idx++)
					gCuts->vals[cnt]->beta[coord->colsC[idx]] += sigma->vals[iStar.sigma].piC[idx];

				for (idx = 1; idx <= num->rvColCnt; idx++)
					gCuts->vals[cnt]->beta[coord->rvCols[idx]] += delta->vals[iStar.delta][observ[obs]].piC[idx];

				count++;
			}
		}

		/* Take the average of the alpha and beta values */
		for (idx = 0; idx <= lenX; idx++)
			gCuts->vals[cnt]->beta[idx] /= (double) k;

		gCuts->vals[cnt]->alpha /= (double) k;

		if (lbType == NONTRIVIAL)
			gCuts->vals[cnt]->alpha += (1 - (double) count / (double) k) * lb;
	}

}//END reform_cuts

/* This function is to calculate the lower bound on the optimal value which is used in stopping rule in full_test() in optimal.c in the case of
 regularized approach. */
double calcBootstrpLB(probType *prob, vectorC incumbX, vectorC piM, vectorC djM, int currIter, double quadScalar, cutsType *cuts) {
	double *bk; 			/* vector: b - A*incumb_x. */
	double *lambda; 		/* vector: the dual of the primal constraints. */
	double bk_lambda; 		/* scalar: bk*lambda. */
	sparseMatrix *A_Trans; 	/* sparse_matrix: the transpose of A(we call it Dbar in our code)*/
	double *A_Trans_lambda; /* vector: - A_Trans * lambda. */
	double theta; 			/* the dual of the reformed cut constraints. */
	double Vk_theta; 		/* Scalar: Vk*theta. */
	double *Bk_theta; 		/* vector: Bk_Transpose * theta, where Bk_Transpose is the matrix of cut coefficients. */
	double *q_vec; 			/* vector: c + Bk_theta - A_Trans_lambda. */
	double q_term; 			/* scalar: q_vec * q_vec. */
	double Lm; 				/* The calculated lower bound of the optimal value. */
	int cnt, i;

	if (!(bk = arr_alloc(prob->num->rows+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to bk", 0);
	if (!(lambda = arr_alloc(prob->num->rows+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to lambda", 0);
	if (!(A_Trans = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans", 0);
	if (!(A_Trans->val = arr_alloc(prob->Dbar->cnt+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->val", 0);
	if (!(A_Trans->row = arr_alloc(prob->Dbar->cnt+1, int)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->row", 0);
	if (!(A_Trans->col = arr_alloc(prob->Dbar->cnt+1, int)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_Trans->col", 0);
	if (!(A_Trans_lambda = arr_alloc(prob->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to A_lambda", 0);
	if (!(Bk_theta = arr_alloc(prob->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to Bk_theta", 0);
	if (!(q_vec = arr_alloc(prob->num->cols+1, double)))
		errMsg("Allocation", "cal_temp_lb", "fail to allocate memory to q_vec", 0);

	/* 1a. Calculate bk, which is A*incumb_x - b. Note: in fact, we are
	 ** calculating -bk here, due to the way function MSparsexvSub works. Also be aware of the one-norm. */
	for (cnt = 0; cnt < prob->num->rows; cnt++)
		bk[cnt + 1] = prob->sp->rhsx[cnt];

	/* 1b. Calculate bk = b - A * incumb_x. */
	MSparsexvSub(prob->Dbar, incumbX, bk);

	/* 1c. Obtain lambda from cell->pi of original master problem constraints. */
	/* Dual values' sign need to be flipped here before assigning to lambda */
	for (cnt = 0; cnt < prob->num->rows; cnt++)
		lambda[cnt + 1] = -piM[cnt + 1];

	/* 1d. Calculate bk_lambda = bk * lambda. */
	bk_lambda = vXv(bk, lambda, NULL, prob->num->rows);

	/* 2a. Calculate A_Trans */
	A_Trans->cnt = prob->Dbar->cnt;

	for (cnt = 1; cnt <= A_Trans->cnt; cnt++) {
		A_Trans->val[cnt] = prob->Dbar->val[cnt];
		A_Trans->row[cnt] = prob->Dbar->col[cnt];
		A_Trans->col[cnt] = prob->Dbar->row[cnt];
	}

	/* 2b. Calculate - A_Trans * lambda. */
	MSparsexvSub(A_Trans, lambda, A_Trans_lambda);

	/* 2c. Calculate -A_trans * lambda - c */
	for (i = 0; i < prob->num->cols; i++)
		A_Trans_lambda[i + 1] += djM[i + 1];

	Vk_theta = 0.0;
	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		/* 3a. Obtain theta from c->pi */
		theta = ((double) (currIter - 1) / (double) cuts->vals[cnt]->cutObs) * piM[cuts->vals[cnt]->rowNum + 1];

		/* 3a. Obtain theta from c->pi */
		Vk_theta += theta*(cuts->vals[cnt]->alpha - vXv(cuts->vals[cnt]->beta, incumbX, NULL, prob->num->cols));

		/* 3b. Calculate Bk_theta = theta * beta */
		for (i = 1; i <= prob->num->cols; i++)
			Bk_theta[i] += theta * cuts->vals[cnt]->beta[i];
	}

	/* 4. Calculate the quadratic vector q_vec = Bk_theta[i] - A_Trans_lambda[i].*/
	for (i = 1; i <= prob->num->cols; i++)
		q_vec[i] = prob->dBar->val[i] - Bk_theta[i] - A_Trans_lambda[i];

	q_term = vXv(q_vec, q_vec, NULL, prob->num->cols);

	/* 5. Calculate the lower bound*/
	Lm = Vk_theta + bk_lambda - q_term / quadScalar / 2.0;

	mem_free(bk);
	mem_free(lambda);
	mem_free(A_Trans->col);
	mem_free(A_Trans->row);
	mem_free(A_Trans->val);
	mem_free(A_Trans);
	mem_free(A_Trans_lambda);
	mem_free(Bk_theta);
	mem_free(q_vec);

	return Lm;
}//END calcBootstrpLB

