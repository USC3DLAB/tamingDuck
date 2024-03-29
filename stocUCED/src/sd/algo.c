/*
 * algo.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */

#include "twoSD.h"
#ifdef __unix__
#include "time.h"
#elif defined(_WIN32) || defined(WIN32)
#include "sys/time.h"
#endif

extern stringC outputDir;
extern configType config;

FILE *file;

int algo(oneProblem *orig, timeType *tim, stocType *stoc, stringC probName, vectorC *edSols, double *objVal) {
	cellType *cell;
	vectorC	 xk = NULL, lb = NULL;
	probType **prob = NULL;
	double	 totalTime;
	FILE 	*soln;
	//struct timeval time;
	
	/* open solver */
	openSolver();
	
	/* open output */
	file = openFile(outputDir, "sd.log", "a");

	/* complete necessary initialization for the algorithm */
	if ( setupAlgo(orig, stoc, tim, &prob, &cell) )
		goto TERMINATE;

	//gettimeofday(&time,NULL);	// gets wall-time
	totalTime = 0; // (double)time.tv_sec + (double)time.tv_usec * .000001;

	/* Use two-stage algorithm to solve the problem */
	if ( solveCell(stoc, prob, cell, probName) ) {
		errMsg("algorithm", "algo", "failed to solve the cells using 2SD algorithm", 0);
		goto TERMINATE;
	}

	//gettimeofday(&time,NULL);
	totalTime = 0; // SMH (double)time.tv_sec + (double)time.tv_usec * .000001 - totalTime;
	
	/* Write solution statistics for optimization process */
	fprintf(file, "\n\nLower bound estimate                   : %f\n", cell->incumbEst);
	fprintf(file, "Total time = %f\n", totalTime);
	//writeStatistic(&soln, prob[0], cell, probName);

	/* evaluating the optimal solution*/
	if (config.EVAL_FLAG == 1) {
		evaluate(&soln, stoc, prob, cell, cell->incumbX);
	}

	fprintf(file, "\n\t\tSuccessfully completed two-stage stochastic decomposition algorithm.\n");

	/* Extract relevant part of the solutions */
	(*edSols) = (vectorC) arr_alloc(prob[0]->num->cols, double);
//	printf("%d variables", prob[0]->num->cols);
	for (int n = 0; n < prob[0]->num->cols; n++ ) {
		(*edSols)[n] = cell->incumbX[n+1];
//		printf(cell->master->cname[n]);
//		printf(" ");
//		printf("%.2f", cell->incumbX[n+1]);
//		printf("\n");
	}
	*objVal = cell->incumbEst;

	/* free up memory before leaving */
	fclose(file);
	if (xk) mem_free(xk);
	if (lb) mem_free(lb);
	if(cell) freeCellType(cell);
	freeProbType(prob, 2);
	return 0;

	TERMINATE:
	fclose(file);
	if(xk) mem_free(xk);
	if(lb) mem_free(lb);
	if(cell) freeCellType(cell);
	if(prob) freeProbType(prob, 2);
	return 1;
}//END algo()

int solveCell(stocType *stoc, probType **prob, cellType *cell, stringC probName) {
	vectorC 	observ;
	int		m, omegaIdx, candidCut;
	CBOOL 	newOmegaFlag;

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Main Algorithm -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	if ( !(observ = (vectorC) arr_alloc(stoc->numOmega + 1, double)) )
		errMsg("allocation", "solveMASP", "observ", 0);

	/******* 0. Initialization: The algorithm begins by solving the master problem as a QP *******/
	while (cell->optFlag == CFALSE && cell->k < config.MAX_ITER) {
		cell->k++;

#if defined(STOCH_CHECK) || defined(ALGO_CHECK)
		printf("\n\t\tIteration-%d :: \n", cell->k);
#else
		if ( (cell->k -1) % 100 == 0) {
			fprintf(file, "\n\t\tIteration-%4d: ", cell->k);
			fflush(file);
		}		
#endif

		/******* 1. Optimality tests *******/
		if (optimal(prob, cell))
			break;

		/******* 2. Generate new observation, and add it to the set of observations *******/
		/* (a) Use the stoc file to generate observations */
		generateOmega(stoc, observ, &config.RUN_SEED);

		/* (b) Since the problem already has the mean values on the right-hand side, remove it from the original observation */
		for ( m = 0; m < stoc->numOmega; m++ ) {
			observ[m] -= stoc->mean[m];
		}

		/* (d) update omegaType with the latest observation. If solving with incumbent then this update has already been processed. */
		omegaIdx = calcOmega(observ - 1, 0, prob[1]->num->numRV, cell->omega, &newOmegaFlag);

		/******* 3. Solve the subproblem with candidate solution, form and update the candidate cut *******/
		if ( (candidCut = formSDCut(prob[1], cell, cell->candidX, omegaIdx, CTRUE)) < 0 ) {
			errMsg("algorithm", "solveCell", "failed to add candidate cut", 0);
			return 1;
		}
		/******* 4. Solve subproblem with incumbent solution, and form an incumbent cut *******/
		if (((cell->k - cell->iCutUpdt) % config.TAU == 0 ) ) {
			if ( (cell->iCutIdx = formSDCut(prob[1], cell, cell->incumbX, omegaIdx, CFALSE) ) < 0 ) {
				errMsg("algorithm", "solveCell", "failed to create the incumbent cut", 0);
				return 1;
			}
			cell->iCutUpdt = cell->k;
		}
		/******* 5. Check improvement in predicted values at candidate solution *******/
		if ( !(cell->incumbChg) && cell->k > 1)
			/* If the incumbent has not changed in the current iteration */
			checkImprovement(prob[0], cell, candidCut);

		/******* 6. Solve the master problem to obtain the new candidate solution */
		if ( solveQPMaster(prob[0]->num, prob[0]->dBar, cell, prob[0]->sp->mar, prob[0]->lb) ) {
			errMsg("algorithm", "solveMASP", "failed to solve master problem", 0);
			return 1;
		}
	}//END while loop

	mem_free(observ);
	return 0;
}//END solveCell()

void writeStatistic(FILE **soln, probType *prob, cellType *cell, stringC probName) {

	(*soln) = openFile(outputDir, "summary.dat", "w");

	fprintf((*soln), "----------------------------------- Problem Information ------------------------------------\n\n");
	fprintf((*soln), "Problem                                : %s\n", probName);
	fprintf((*soln), "First Stage Rows                       : %d\n", prob->num->rows);
	fprintf((*soln), "First Stage Columns                    : %d\n", prob->num->cols);

	fprintf((*soln), "\n--------------------------------------- Optimization ---------------------------------------\n\n");

	fprintf((*soln), "Algorithm                              : Two-stage Stochastic Decomposition\n");
	fprintf((*soln), "Number of iterations                   : %d\n", cell->k);
	fprintf((*soln), "Lower bound estimate                   : %f\n", cell->incumbEst);
	fclose((*soln));

}//END WriteStat
