/*
 * prob.c
 *
 *  Created on: Sep 30, 2015
 *      Author: Harsha Gangammanavar
 */

#include "../solverUtilities/utils.h"
#include "../solverUtilities/solver.h"
#include "smps.h"
#include "prob.h"

/* Decomposes the problem _orig_ into subproblems as well as decomposes the stochastic information _stoc_ into stage stochastic information. The decomposition
 * is carried out using information specified in _tim_. The function also stores stage lower bound information provided in _Lb_. It returns an array of
 * probType structures, each probType corresponds to a particular stage */
probType **newProb(oneProblem *orig, stocType *stoc, timeType *tim, vectorC lb, double TOLERANCE) {
	probType **prob;
	char	 *q;
	int		 i, k, m, t, rOffset = 0, cOffset = 0;

	/* allocate memory to elements of probType */
	if ( !(prob = (probType **) arr_alloc(tim->numStages, probType *)) )
		errMsg("allocation", "newProb", "prob", 0);

	/* allocate memory to members of probType for stagewise subProblems, and allocate values to static fields*/
	for ( t = 0; t < tim->numStages; t++ ) {
		if ( !(prob[t] = (probType *) mem_malloc(sizeof(probType))) )
			errMsg("allocation", "newProb", "prob[t]", 0);
		if ( !(prob[t]->name = (stringC) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "newProb", "stage name", 0);
		if( !(prob[t]->sp = (oneProblem *) mem_malloc (sizeof(oneProblem))))
			errMsg("allocation", "newProb", "stage problem", 0);
		prob[t]->sp->lp = NULL; prob[t]->lb = 0.0;

		if (  !(prob[t]->dBar = (sparseVector *) mem_malloc(sizeof(sparseVector))) )
			errMsg("allocation", "newProb", "stage problem cost coefficients", 0);
		if ( !(prob[t]->dBar->col = (intvec) arr_alloc(orig->mac+1, int)) )
			errMsg("allocation", "newProb", "stage problem cost coefficients columns", 0);
		if ( !(prob[t]->dBar->val = (vectorC) arr_alloc(orig->mac+1, double)) )
			errMsg("allocation", "newProb", "stage problem cost coefficients values", 0);
		prob[t]->dBar->cnt = 0;

		if ( !(prob[t]->bBar = (sparseVector *) mem_malloc(sizeof(sparseVector))) )
			errMsg("allocation", "newProb", "stage problem right hand side", 0);
		if ( !(prob[t]->bBar->col = (intvec) arr_alloc(orig->mar+1, int)) )
			errMsg("allocation", "newProb", "stage problem right-hand rows", 0);
		if ( !(prob[t]->bBar->val = (vectorC) arr_alloc(orig->mar+1, double)) )
			errMsg("allocation", "newProb", "stage problem right-hand values", 0);
		prob[t]->bBar->cnt = 0;

		strcpy(prob[t]->name, tim->stgNames[t]);

		if ( t < tim->numStages - 1 ) {
			prob[t]->sp->mar = prob[t]->sp->marsz = tim->row[t+1] - tim->row[t];
			prob[t]->sp->mac = prob[t]->sp->macsz = tim->col[t+1] - tim->col[t];
			prob[t]->sp->rstorsz = orig->rname[tim->row[t+1]] - orig->rname[tim->row[t]];
			prob[t]->sp->cstorsz = orig->cname[tim->col[t+1]] - orig->cname[tim->col[t]];
			rOffset += prob[t]->sp->rstorsz;
			cOffset += prob[t]->sp->cstorsz;
		}
		else {
			prob[t]->sp->mar = prob[t]->sp->marsz = orig->mar - tim->row[t];
			prob[t]->sp->mac = prob[t]->sp->macsz = orig->mac - tim->col[t];
			prob[t]->sp->rstorsz = orig->rstorsz - rOffset;
			prob[t]->sp->cstorsz = orig->cstorsz - cOffset;
		}
		prob[t]->sp->numInt = 0;
		prob[t]->sp->numBin = 0;
		prob[t]->sp->matsz = 0;
		prob[t]->sp->numnz = 0;
		prob[t]->sp->objsen = orig->objsen;
		prob[t]->sp->type = PROB_LP;

		/* stage oneProblem */
		if(!(prob[t]->sp->name = (stringC) arr_alloc(NAMESIZE, char)))
			errMsg("allocation", "newProb", "stage problem name", 0);
		if(!(prob[t]->sp->objname = (stringC) arr_alloc(NAMESIZE, char)))
			errMsg("allocation", "newProb", "stage problem objname", 0);
		if(!(prob[t]->sp->objx = (vectorC) arr_alloc(prob[t]->sp->macsz, double)))
			errMsg("allocation", "newProb", "stage problem objx", 0);
		if(!(prob[t]->sp->bdl = (vectorC) arr_alloc(prob[t]->sp->macsz, double)))
			errMsg("allocation", "newProb", "stage problem bdl", 0);
		if(!(prob[t]->sp->bdu = (vectorC) arr_alloc(prob[t]->sp->macsz, double)))
			errMsg("allocation", "newProb", "stage problem bdu", 0);
		if(!(prob[t]->sp->ctype = (stringC) arr_alloc(prob[t]->sp->macsz, char)))
			errMsg("allocation", "newProb", "stage problem column type", 0);
		if(!(prob[t]->sp->rhsx = (vectorC) arr_alloc(prob[t]->sp->marsz, double)))
			errMsg("allocation", "newProb", "stage problem rhsx", 0);
		if(!(prob[t]->sp->senx = (stringC) arr_alloc(prob[t]->sp->marsz, char)))
			errMsg("allocation", "newProb", "stage problem senx", 0);
		if(!(prob[t]->sp->matbeg = (intvec) arr_alloc(prob[t]->sp->macsz, int)))
			errMsg("allocation", "newProb", "stage problem matbeg", 0);
		if(!(prob[t]->sp->matcnt = (intvec) arr_alloc(prob[t]->sp->macsz, int)))
			errMsg("allocation", "newProb", "stage problem matcnt", 0);
		if(!(prob[t]->sp->cname = (stringC *) arr_alloc(prob[t]->sp->macsz, stringC)))
			errMsg("allocation", "newProb", "stage problem cname", 0);
		if(!(prob[t]->sp->cstore = (stringC) arr_alloc(prob[t]->sp->cstorsz, char)))
			errMsg("allocation", "newProb", "stage problem cstore", 0);
		if(!(prob[t]->sp->rname = (stringC *) arr_alloc(prob[t]->sp->marsz, stringC)))
			errMsg("allocation", "newProb", "stage problem rname", 0);
		if(!(prob[t]->sp->rstore = (stringC) arr_alloc(prob[t]->sp->rstorsz, char)))
			errMsg("allocation", "newProb", "stage problem rstore", 0);
		if(!(prob[t]->sp->matval = (vectorC) arr_alloc(orig->matsz, double)))
			errMsg("allocation", "newProb", "stage problem matval", 0);
		if(!(prob[t]->sp->matind = (intvec) arr_alloc(orig->matsz, int)))
			errMsg("allocation", "newProb", "stage problem matind", 0);
		strcpy(prob[t]->sp->objname, orig->objname);
		sprintf(prob[t]->sp->name, "%s_%d", orig->name, t);

		/* stage transfer matrix */
		if ( t == 0 )
			prob[t]->Cbar = NULL;
		else {
			if ( !(prob[t]->Cbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
				errMsg("allocation", "newProb", "stage transfer matrix", 0);
			if(!(prob[t]->Cbar->row = (intvec) arr_alloc(orig->matsz + 1, int)))
				errMsg("allocation", "newProb", "transfer matrix rows", 0);
			if(!(prob[t]->Cbar->col = (intvec) arr_alloc(orig->matsz + 1, int)))
				errMsg("allocation", "newProb", "transfer matrix columns", 0);
			if(!(prob[t]->Cbar->val = (vectorC) arr_alloc(orig->matsz + 1, double)))
				errMsg("allocation", "newProb", "transfer matrix values", 0);
			prob[t]->Cbar->cnt = 0;
		}

		/* TODO (HG): stage dynamics: include dynamics to the model */
		if ( t == 0) {
			prob[t]->Abar = NULL;
			prob[t]->Bbar = NULL;
			prob[t]->aBar = NULL;
			prob[t]->cBar = NULL;
		}
		else {
			prob[t]->Abar = NULL;
			prob[t]->Bbar = NULL;
			prob[t]->aBar = NULL;
			prob[t]->cBar = NULL;
		}

		/* Stage recourse matrix: Dbar is used to setup the QP in the forward pass. Since the QP is used only for non-terminal
		 * stages, Dbar for terminal stage is set to NULL */
		if ( t == tim->numStages-1 )
			prob[t]->Dbar = NULL;
		else {
			if ( !(prob[t]->Dbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
				errMsg("allocation", "newProb", "stage constraint matrix", 0);
			if(!(prob[t]->Dbar->row = (intvec) arr_alloc(orig->matsz+1, int)))
				errMsg("allocation", "newProb", "stage constraint matrix rows", 0);
			if(!(prob[t]->Dbar->col = (intvec) arr_alloc(orig->matsz+1, int)))
				errMsg("allocation", "newProb", "stage constraint matrix columns", 0);
			if(!(prob[t]->Dbar->val = (vectorC) arr_alloc(orig->matsz+1, double)))
				errMsg("allocation", "newProb", "stage constraint matrix values", 0);
			prob[t]->Dbar->cnt = 0;
		}
	}

	for ( t = 0; t < tim->numStages-1; t++ ) {
		/* lower bound on cost-to-go function */
		if ( lb != NULL )
			prob[t]->lb = lb[t];

		/* copy column names of non-terminal stage*/
		m = 0;
		for ( q = orig->cname[tim->col[t]]; q < orig->cname[tim->col[t+1]]; q++ )
			prob[t]->sp->cstore[m++] = *q;

		/* copy column information for non-terminal stage */
		cOffset = prob[t]->sp->cstore - orig->cname[tim->col[t]];
		for ( m = tim->col[t]; m < tim->col[t+1]; m++ ) {
			k = m - tim->col[t];
			prob[t]->dBar->val[prob[t]->dBar->cnt+1] = orig->objx[m];
			prob[t]->dBar->col[prob[t]->dBar->cnt+1] = m - tim->col[t]+1;
			prob[t]->dBar->cnt++;
			prob[t]->sp->objx[k] = orig->objx[m];
			prob[t]->sp->bdl[k] = orig->bdl[m];
			prob[t]->sp->bdu[k] = orig->bdu[m];
			if ( orig->ctype[m] == 'I')
				prob[t]->sp->numInt++;
			else if ( orig->ctype[m] == 'B' )
				prob[t]->sp->numBin++;
			prob[t]->sp->ctype[k] = orig->ctype[m];
			prob[t]->sp->cname[k] = orig->cname[m] + cOffset;
			prob[t]->sp->matcnt[k] = 0;
			for ( i = orig->matbeg[m]; i < orig->matbeg[m]+orig->matcnt[m]; i++ ) {
				if (orig->matind[i] < tim->row[t+1]) {
					/* The coefficient is part of the current stage constraint matrix */
					if ( k == 0 )
						prob[t]->sp->matbeg[k] = 0;
					else
						prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1] + prob[t]->sp->matcnt[k-1];
					prob[t]->sp->matval[prob[t]->sp->matsz] = orig->matval[i];
					prob[t]->sp->matind[prob[t]->sp->matsz] = orig->matind[i] - tim->row[t];
					++prob[t]->sp->matcnt[k];
					++prob[t]->sp->matsz;
					++prob[t]->sp->numnz;
					prob[t]->Dbar->val[prob[t]->Dbar->cnt+1] = orig->matval[i];
					prob[t]->Dbar->col[prob[t]->Dbar->cnt+1] = m-tim->col[t]+1;
					prob[t]->Dbar->row[prob[t]->Dbar->cnt+1] = orig->matind[i] - tim->row[t]+1;
					++prob[t]->Dbar->cnt;
				}
				else {
					/* The coefficient is part of the next stage transfer matrix */
					prob[t+1]->Cbar->val[prob[t+1]->Cbar->cnt+1] = orig->matval[i];
					prob[t+1]->Cbar->row[prob[t+1]->Cbar->cnt+1] = orig->matind[i] - tim->row[t+1]+1;
					prob[t+1]->Cbar->col[prob[t+1]->Cbar->cnt+1] = m-tim->col[t]+1;
					++prob[t+1]->Cbar->cnt;
				}
			}
		}

		/* if integer or binary variables are encountered, then label the stage problem as a mixed integer LP */
		if ( prob[t]->sp->numInt + prob[t]->sp->numBin > 0 )
			prob[t]->sp->type = PROB_MILP;

		/* copy row name of non-terminal stage */
		m = 0;
		for ( q = orig->rname[tim->row[t]]; q < orig->rname[tim->row[t+1]]; q++ )
			prob[t]->sp->rstore[m++] = *q;

		/* copy row information for non-terminal stage */
		rOffset = prob[t]->sp->rstore - orig->rname[tim->row[t]];
		for ( m = tim->row[t]; m < tim->row[t+1]; m++ ) {
			k = m - tim->row[t];
			prob[t]->sp->rhsx[k] = orig->rhsx[m];
			prob[t]->sp->senx[k] = orig->senx[m];
			prob[t]->sp->rname[k] = orig->rname[m]+rOffset;
			prob[t]->bBar->val[prob[t]->bBar->cnt+1] = orig->rhsx[m];
			prob[t]->bBar->col[prob[t]->bBar->cnt+1] = m - tim->row[t]+1;
			prob[t]->bBar->cnt++;
		}
	}

	/* Now copy the terminal stage problem */
	/* copy column information for terminal stage*/
	m = 0;
	for ( q = orig->cname[tim->col[t]]; q < orig->cname[0] + orig->cstorsz; q++ )
		prob[t]->sp->cstore[m++] = *q;

	cOffset = prob[t]->sp->cstore - orig->cname[tim->col[t]];
	for ( m = tim->col[t]; m < orig->mac; m++ ) {
		k = m - tim->col[t];
		prob[t]->dBar->val[prob[t]->dBar->cnt+1] = orig->objx[m];
		prob[t]->dBar->col[prob[t]->dBar->cnt+1] = m-tim->col[t]+1;
		prob[t]->dBar->cnt++;
		prob[t]->sp->objx[k] = orig->objx[m];
		prob[t]->sp->bdl[k] = orig->bdl[m];
		prob[t]->sp->bdu[k] = orig->bdu[m];
		if (orig->ctype[m] != 'C') {
			errMsg("setup", "newProb", "integer variable in non-root stage",0);
			return NULL;
		}
		else
			prob[t]->sp->ctype[k] = orig->ctype[m];
		prob[t]->sp->cname[k] = orig->cname[m] + cOffset;
		prob[t]->sp->matcnt[k] = 0;
		if ( orig->matcnt[m] > 0 )
			for ( i = orig->matbeg[m]; i < orig->matbeg[m]+orig->matcnt[m]; i++ ) {
				/* The coefficient is part of the current stage constraint matrix */
				if ( k == 0 )
					prob[t]->sp->matbeg[k] = 0;
				else
					prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1] + prob[t]->sp->matcnt[k-1];
				prob[t]->sp->matval[prob[t]->sp->matsz] = orig->matval[i];
				prob[t]->sp->matind[prob[t]->sp->matsz] = orig->matind[i]-tim->row[t];
				++prob[t]->sp->matcnt[k];
				++prob[t]->sp->matsz;
				++prob[t]->sp->numnz;
			}
		else {
			if ( k == 0 )
				prob[t]->sp->matbeg[k] = 0;
			else
				prob[t]->sp->matbeg[k] = prob[t]->sp->matbeg[k-1];
		}
	}

	/* copy row information for terminal stage */
	m = 0;
	for ( q = orig->rname[tim->row[t]]; q < orig->rname[0] + orig->rstorsz; q++ )
		prob[t]->sp->rstore[m++] = *q;

	rOffset = prob[t]->sp->rstore - orig->rname[tim->row[t]];
	for ( m = tim->row[t]; m < orig->mar; m++ ) {
		k = m - tim->row[t];
		prob[t]->sp->rhsx[k] = orig->rhsx[m];
		prob[t]->sp->senx[k] = orig->senx[m];
		prob[t]->sp->rname[k] = orig->rname[m]+rOffset;
		prob[t]->bBar->val[prob[t]->bBar->cnt+1] = orig->rhsx[m];
		prob[t]->bBar->col[prob[t]->bBar->cnt+1] = m - tim->row[t]+1;
		prob[t]->bBar->cnt++;
	}

#ifdef DECOMPOSE_CHECK
	/* write stage problems in LP format to verify decomposition */
	char fname[BLOCKSIZE];
	for ( t = 0; t < tim->numStages; t++) {
		if ( !(prob[t]->sp->lp = setupProblem(prob[t]->sp->name, prob[t]->sp->type, prob[t]->sp->mac, prob[t]->sp->mar, prob[t]->sp->objsen, prob[t]->sp->objx, prob[t]->sp->rhsx, prob[t]->sp->senx,
				prob[t]->sp->matbeg, prob[t]->sp->matcnt, prob[t]->sp->matind, prob[t]->sp->matval, prob[t]->sp->bdl, prob[t]->sp->bdu, NULL, prob[t]->sp->cname, prob[t]->sp->rname, prob[t]->sp->ctype)) ) {
			errMsg("solver", "newProb", "failed to setup stage problem in solver", 0);
			return prob;
		}
		sprintf(fname, "stageProb%d.lp", t);
		if ( writeProblem(prob[t]->sp->lp, fname) ) {
			errMsg("solver", "newProb", "failed to write stage problem", 0);
			return prob;
		}
	}
#endif

	/* save size information in numType */
	for ( t = 0; t < tim->numStages; t++ ) {
		if ( !(prob[t]->num = (numType *) mem_malloc(sizeof(numType))) )
			errMsg("allocation", "newProb", "prob[t]->num",0);
		prob[t]->num->cols = prob[t]->sp->mac;
		prob[t]->num->rows = prob[t]->sp->mar;
		if ( t == 0 )
			prob[t]->num->prevCols = prob[t]->num->prevRows = 0;
		else {
			prob[t]->num->prevCols = prob[t-1]->num->cols;
			prob[t]->num->prevRows = prob[t-1]->num->rows;
		}
	}

	/* decompose the stochastic elements of the problem */
	for ( t = 0; t < tim->numStages; t++ ) {
		if ( t == 0)
			prob[t]->omegas = NULL;
		else {
			if ( !(prob[t]->omegas = (omegastuff *) mem_malloc(sizeof(omegastuff))) )
				errMsg("allocation", "newProb", "prob->omegas", 0);
			prob[t]->omegas->col = NULL; prob[t]->omegas->mean = NULL; prob[t]->omegas->row = NULL;
			prob[t]->omegas->numRV = 0; prob[t]->omegas->beg = -1;
		}
		prob[t]->num->numRV = prob[t]->num->rvColCnt = prob[t]->num->rvRowCnt = 0;
		prob[t]->num->rvAOmCnt = prob[t]->num->rvBOmCnt = prob[t]->num->rvCOmCnt = prob[t]->num->rvDOmCnt = 0;
		prob[t]->num->rvaOmCnt = prob[t]->num->rvbOmCnt = prob[t]->num->rvcOmCnt = prob[t]->num->rvdOmCnt = 0;
		prob[t]->num->intCols = prob[t]->sp->numInt;
		prob[t]->num->binCols = prob[t]->sp->numBin;
	}

	/* go through the list of random variable and assign them to appropriate parts (right-hand side and objective coefficients) */
	for ( m = 0; m < stoc->numOmega; m++ ) {
		if ( stoc->col[m] == -1 ) {
			/* randomness in right-hand side */
			t = 0;
			while ( t < tim->numStages ) {
				if ( stoc->row[m] < tim->row[t] )
					break;
				t++;
			}
			t--;
		}
		else {
			/* randomness in either objective function coefficients or the transfer matrix */
			t = 0;
			while ( t < tim->numStages ) {
				if ( stoc->col[m] < tim->col[t] )
					break;
				t++;
			}
			t--;
		}

		if ( t == 0 ) {
			errMsg("setup", "newProb", "encountered randomness is root-stage", 0);
			return NULL;
		}

		if ( prob[t]->omegas->numRV == 0) {
			if ( !(prob[t]->omegas->col = (intvec) arr_alloc(stoc->numOmega+1, int)) )
				errMsg("allocation", "newProb", "prob->omegas->col", 0);
			if ( !(prob[t]->omegas->row = (intvec) arr_alloc(stoc->numOmega+1, int)) )
				errMsg("allocation", "newProb", "prob->omegas->row", 0);
			if ( !(prob[t]->omegas->mean = (vectorC) arr_alloc(stoc->numOmega+1, double)) )
				errMsg("allocation", "newProb", "prob->omegas->mean", 0);
			prob[t]->omegas->beg = m;
		}

		if ( stoc->col[m] == -1 )
			prob[t]->omegas->col[prob[t]->omegas->numRV+1] = -1;
		else
			prob[t]->omegas->col[prob[t]->omegas->numRV+1] = stoc->col[m] - tim->col[t]+1;
		if ( stoc->row[m] == -1 )
			prob[t]->omegas->row[prob[t]->omegas->numRV+1] = -1;
		else
			prob[t]->omegas->row[prob[t]->omegas->numRV+1] = stoc->row[m] - tim->row[t]+1;
		prob[t]->omegas->mean[prob[t]->omegas->numRV+1] = stoc->mean[m];
		prob[t]->omegas->numRV++;
		prob[t]->num->numRV++;
		if ( stoc->row[m] == -1 )
			/* random variable in objective coefficient */
			prob[t]->num->rvdOmCnt++;
		if ( stoc->col[m] == -1)
			/* random variable in right-hand side */
			prob[t]->num->rvbOmCnt++;
		else if ( stoc->col[m] < tim->col[t] )
			/* random variable in transfer matrix */
			prob[t]->num->rvCOmCnt++;
	}

	for ( t = 1; t < tim->numStages; t++ ) {
		if (prob[t]->omegas->numRV == 0) {
			freeOmegastuff(prob[t]->omegas);
			prob[t]->omegas = NULL;
		}
		else {
			prob[t]->omegas->col = (intvec) mem_realloc(prob[t]->omegas->col, (prob[t]->omegas->numRV+1)*sizeof(int));
			prob[t]->omegas->row = (intvec) mem_realloc(prob[t]->omegas->row, (prob[t]->omegas->numRV+1)*sizeof(int));
			prob[t]->omegas->mean = (vectorC) mem_realloc(prob[t]->omegas->mean, (prob[t]->omegas->numRV+1)*sizeof(double));
		}
	}

	/* save coordinate information in coordType */
	prob[0]->coord = NULL;
	for ( t = 1; t < tim->numStages; t++ ) {
		if ( !(prob[t]->coord = (coordType *) mem_malloc(sizeof(coordType))) )
			errMsg("allocation", "newProb", "prob[t]->coord",0);

		if ( prob[t]->omegas != NULL ) {
			if ( !(prob[t]->coord->omegaCol = (intvec) arr_alloc(prob[t]->num->numRV+1, int)) )
				errMsg("allocation", "newProb", "prob[t]->coord->omegaCol", 0);
			if ( !(prob[t]->coord->omegaRow = (intvec) arr_alloc(prob[t]->num->numRV+1, int)) )
				errMsg("allocation", "newProb", "prob[t]->coord->omegaCol", 0);
			for ( m = 1; m <= prob[t]->num->numRV; m++ ) {
				prob[t]->coord->omegaCol[m] = prob[t]->omegas->col[m];
				prob[t]->coord->omegaRow[m] = prob[t]->omegas->row[m];
			}
			prob[t]->coord->rvCols = findElems(prob[t]->coord->omegaCol, prob[t]->num->numRV, &prob[t]->num->rvColCnt);
			prob[t]->coord->rvRows = findElems(prob[t]->coord->omegaRow, prob[t]->num->numRV, &prob[t]->num->rvRowCnt);
		}
		else {
			prob[t]->coord->omegaCol = NULL;
			prob[t]->coord->omegaRow = NULL;
			prob[t]->coord->rvCols = NULL;
			prob[t]->coord->rvRows = NULL;
		}
		prob[t]->coord->colsC = findElems(prob[t]->Cbar->col, prob[t]->Cbar->cnt, &prob[t]->num->cntCcols);
	}

	/* modify the bBar and Cbar with mean values computed from stoch file */
	cOffset = 0;
	for ( t = 1; t < tim->numStages; t++ ) {
		for ( m = 1; m <= prob[t]->num->numRV; m++ ) {
			if ( prob[t]->coord->omegaCol[m] == -1) {
				i = 1;
				while ( i <= prob[t]->bBar->cnt ) {
					if ( prob[t]->bBar->col[i] == prob[t]->coord->omegaRow[m])
						break;
					i++;
				}
				prob[t]->bBar->val[i] = stoc->mean[cOffset+m-1];
			}
			else {
				i = 1;
				while ( i <= prob[t]->Cbar->cnt ) {
					if ( prob[t]->Cbar->col[i] == prob[t]->coord->omegaCol[m] && prob[t]->Cbar->row[i] == prob[t]->coord->omegaRow[m] )
						break;
					i++;
				}
				prob[t]->Cbar->val[i] = stoc->mean[cOffset+m-1];
			}
		}
		cOffset += prob[t]->num->numRV;
	}

#if 1
	t = 1;
	prob[t]->sp->lp = setupProblem(prob[t]->sp->name, prob[t]->sp->type, prob[t]->sp->mac, prob[t]->sp->mar, prob[t]->sp->objsen, prob[t]->sp->objx, prob[t]->sp->rhsx, prob[t]->sp->senx,
			prob[t]->sp->matbeg, prob[t]->sp->matcnt, prob[t]->sp->matind, prob[t]->sp->matval, prob[t]->sp->bdl, prob[t]->sp->bdu, NULL, prob[t]->sp->cname, prob[t]->sp->rname,
			prob[t]->sp->ctype);
	if ( prob[t]->sp->lp == NULL ) {
		errMsg("solver", "newSubprob", "subprob",0);
		return NULL;
	}
#endif

	return prob;
}//END newProb()

/* setup and solve the original problem _orig_ with expected values for all random variables provided in _stoc_. If the problem is an mixed-integer program,
 *  then a relaxed problem is solved. The function returns a vector of mean value solutions, if there is an error it returns NULL.*/
vectorC meanProblem(oneProblem *orig, stocType *stoc) {
	vectorC	xk;
	double	obj = 0.0;
	int 	n, status;

	/* setup problem in the solver */
	orig->lp = setupProblem(orig->name, orig->type, orig->mac, orig->mar, orig->objsen, orig->objx, orig->rhsx, orig->senx, orig->matbeg, orig->matcnt,
			orig->matind, orig->matval, orig->bdl, orig->bdu, NULL, orig->cname, orig->rname, orig->ctype);
	if ( orig->lp == NULL ) {
		errMsg("setup", "meanProblem", "failed to setup the mean problem", 0);
		return NULL;
	}

	/* change the coefficients and right-hand side to mean values */
	for (n = 0; n < stoc->numOmega; n++ ) {
		status = changeCoef(orig->lp, stoc->row[n], stoc->col[n], stoc->mean[n]);
		if ( status ) {
			errMsg("setup", "meanProblem", "failed to change the coefficients with mean values", 0);
			return NULL;
		}
	}

	/* change the problem type to solve the relaxed mean value problem */
	if ( orig->type == PROB_MILP || orig->type == PROB_MIQP ) {
		status = changeProbType(orig->lp, PROB_LP);
		if ( status ) {
			errMsg("solver", "meanProblem", "failed to relax the mixed-integer program", 0);
			return NULL;
		}
	}

	/* write the mean value problem */
//	status = writeProblem(orig->lp, "original.lp");
//	if ( status ) {
//		errMsg("solver", "meanProblem", "failed to write the problem", 0);
//		return NULL;
//	}

	/* solve the mean value problem */
	//MARK: PROB_LP is changed to PROB_Qp
	status = solveProblem(orig->lp, orig->name, PROB_QP, &status);
	if ( status ) {
		errMsg("setup", "meanProblem", "failed to solve mean value problem", 0);
		return NULL;
	}

	/* obtain solution information and print */
	if ( !(xk = (vectorC) arr_alloc(orig->mac+1, double)) )
		errMsg("allocation", "meanProblem", "sol", 0);

	/* print results */
	obj = getObjective(orig->lp, PROB_LP);
#if defined (VERBOSE)
	printf("Optimal objective function value for (relaxed) mean value problem = %lf\n", obj);
#endif

	/* obtain the primal solution */
	getPrimal(orig->lp,	xk, orig->mac);

	return xk;
}//END meanProblem()

vectorC calcLowerBound(oneProblem *orig, timeType *tim, stocType *stoc) {
	sparseVector	*bBar;
	sparseMatrix	*Cbar;
	vectorC		duals, vals, beta, lb;
	intvec		indices;
	double		alpha;
	int 		status, stat1, t, col, row, m, n;
	LPptr		lpClone;
	BOOL		zeroLB;

	if ( !(lb = (vectorC) arr_alloc(tim->numStages, double)) )
		errMsg("allocation", "getLowerBound", "lb", 0);
	if ( !(duals = (vectorC) arr_alloc(orig->mar+1, double)) )
		errMsg("allocation", "getLowerBound", "duals", 0);
	if (!(indices = (intvec) arr_alloc(orig->mac, int)))
		errMsg("allocation", "getLowerBound", "indices", 0);
	if (!(vals = (vectorC) arr_alloc(orig->mac, double)))
		errMsg("allocation", "getLowerBound", "vals", 0);
	if (!(bBar = (sparseVector *) mem_malloc(sizeof(sparseVector))) )
		errMsg("allocation", "getLowerBound", "bBar", 0);
	if (!(bBar->col = (intvec) arr_alloc(orig->mar+1, int)) )
		errMsg("allocation", "getLowerBound", "bBar->col", 0);
	if (!(bBar->val = (vectorC) arr_alloc(orig->mar+1, double)) )
		errMsg("allocation", "getLowerBound", "bBar->val", 0);
	if (!(Cbar = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
		errMsg("allocation", "getLowerBound", "Cbar", 0);
	if (!(Cbar->col = (intvec) arr_alloc(orig->matsz+1, int)) )
		errMsg("allocation", "getLowerBound", "Cbar->col", 0);
	if (!(Cbar->row = (intvec) arr_alloc(orig->matsz+1, int)) )
		errMsg("allocation", "getLowerBound", "Cbar->row", 0);
	if (!(Cbar->val = (vectorC) arr_alloc(orig->matsz+1, double)) )
		errMsg("allocation", "getLowerBound", "Cbar->val", 0);
	if ( !(beta = (vectorC) arr_alloc(orig->mac+3, double)) )
		errMsg("allocation", "getLowerBound", "beta", 0);

	/* obtain dual solutions from the mean value solve */
	status = getDual(orig->lp, duals, orig->mar);
	if ( status ) {
		errMsg("setup", "getLowerBound", "failed to obtain dual for the mean value problem", 0);
		return NULL;
	}

#if defined (VERBOSE)
	printf("Lower bounds computed = ");
#endif

	for ( t = 1; t < tim->numStages; t++ ) {
		zeroLB = CTRUE;
		col = tim->col[t]; row = tim->row[t];
		bBar->cnt = 0; Cbar->cnt = 0;

		/* check to see if zero is a valid lower bound. If so, then we will use it as a lower bound */
		n = col;
		while ( n < orig->mac ) {
			if ( orig->objx[n]*orig->bdl[n] < 0 || orig->objx[n]*orig->bdu[n] < 0) {
				zeroLB = CFALSE;
				break;
			}
			n++;
		}
		if ( !(zeroLB) ) {
			/* clone to problem to be used for computing the lower bound */
			lpClone = cloneProblem(orig->lp);

			/* extract bBar */
			m = 0;
			for (n = row; n < orig->mar; n++) {
				/* if the element has randomness in right-hand side, then make sure it is accounted */
				if ( n == stoc->row[m] && m < stoc->numOmega)
					bBar->val[bBar->cnt + 1] = stoc->mean[m++];
				else
					bBar->val[bBar->cnt + 1] = orig->rhsx[n];
				bBar->col[bBar->cnt + 1] = n - row + 1;
				++bBar->cnt;
			}

			/* extract Cbar */
			for (n = 0; n < col; n++) {
				for (m = orig->matbeg[n]; m < orig->matbeg[n] + orig->matcnt[n]; m++) {
					if (orig->matind[m] >= row) {
						/* The coefficient is part of the subproblem's T matrix */
						Cbar->val[Cbar->cnt + 1] = orig->matval[m];
						Cbar->row[Cbar->cnt + 1] = orig->matind[m] - row + 1;
						Cbar->col[Cbar->cnt + 1] = n + 1;
						++Cbar->cnt;
					}
				}
			}

			/* compute alpha and beta */
			alpha = 0.0;
			for ( n = 1; n <= bBar->cnt; n++ )
				alpha += duals[bBar->col[n]+row] * bBar->val[n];
			for (n = 0; n <= col + 2; n++)
				beta[n] = 0.0;
			for (n = 1; n <= Cbar->cnt; n++)
				beta[Cbar->col[n]] = beta[Cbar->col[n]] + duals[row + Cbar->row[n]] * Cbar->val[n];
			beta[col + 2] = -1;

			/* use mean-cut coefficients to create the lower-bounding problem */
			for (n = 0; n < col; n++) {
				indices[n] = n;
				vals[n] = -beta[n + 1];
			}
			for (n = col; n < orig->mac; n++) {
				indices[n] = n;
				vals[n] = 0.0;
			}

			/* Change the objective in solver to prepare for lower bound calculation */
			status = changeObjx(lpClone, orig->mac, indices, vals);
			if ( status ) {
				errMsg("setup", "calcLowerBound", "failed to change objective coefficients in solver while computing lower bound", 0);
				return NULL;
			}

#ifdef SETUP_CHECK
			writeProblem(lpClone, "lowerBoundCalc.lp");
#endif
			/* solve the problem */
			status = solveProblem(lpClone, "lowerBoundCalc", PROB_LP, &stat1);
			if ( status ) {
				errMsg("setup", "calcLowerBound", "failed to solve problem computing lower bound", 0);
				return NULL;
			}

			/* get lower bound */
			lb[t-1] = getObjective(lpClone, PROB_LP) + alpha;

			/* release the problem */
			status = freeProblem(lpClone);
			if ( status ) {
				errMsg("setup", "calcLowerBound", "failed to free problem", 0);
				return NULL;
			}
		}
		else
			lb[t-1] = 0.0;

#if defined (VERBOSE)
		printf("%0.3lf\t", lb[t-1]);
#endif
	}
	lb[t-1] = 0.0;

#if defined (VERBOSE)
	printf("\n");
#endif
	
	mem_free(indices);
	mem_free(vals);
	mem_free(beta);
	mem_free(duals);
	freeSparseVector(bBar); freeSparseMatrix(Cbar);

	/* change the problem type back to its original form */
	if ( orig->type == PROB_MILP ) {
		status = changeProbType(orig->lp, PROB_MILP);
		if ( status ) {
			errMsg("solver", "calcLowerBound", "failed to relax the mean value problem", 0);
			return NULL;
		}
	}

#if 0
	for ( t = 1; t < tim->numStages; t++ )
		lb[t] = lb[0];
#endif

	return lb;
}//END calcLowerBound()

/* free up the probType. Needs number of stages as input */
void freeProbType(probType **prob, int T) {
	int t;

	if ( prob ) {
		for ( t = 0; t < T; t++ ) {
			if (prob[t]) {
				if ( prob[t]->Abar ) freeSparseMatrix(prob[t]->Abar);
				if ( prob[t]->Bbar ) freeSparseMatrix(prob[t]->Bbar);
				if ( prob[t]->Cbar ) freeSparseMatrix(prob[t]->Cbar);
				if ( prob[t]->Dbar ) freeSparseMatrix(prob[t]->Dbar);
				if ( prob[t]->aBar ) freeSparseVector(prob[t]->aBar);
				if ( prob[t]->bBar ) freeSparseVector(prob[t]->bBar);
				if ( prob[t]->cBar ) freeSparseVector(prob[t]->cBar);
				if ( prob[t]->dBar ) freeSparseVector(prob[t]->dBar);
				if ( prob[t]->sp ) freeOneProblem(prob[t]->sp);
				if ( prob[t]->name) mem_free(prob[t]->name);
				if ( prob[t]->num ) mem_free(prob[t]->num);
				if ( prob[t]->coord) freeCoordType(prob[t]->coord);
				if ( prob[t]->omegas) freeOmegastuff(prob[t]->omegas);
				mem_free(prob[t]);
			}
		}
		mem_free(prob);
	}

}//END freeProb()

/* free up the coordType */
void freeCoordType (coordType *coord) {

	if (coord->omegaCol) mem_free(coord->omegaCol);
	if (coord->omegaRow) mem_free(coord->omegaRow);
	if (coord->colsC) mem_free(coord->colsC);
	if (coord->rvCols) mem_free(coord->rvCols);
	if (coord->rvRows) mem_free(coord->rvRows);

	mem_free(coord);

}//END freeCoordType()

/* free up omegaStuff */
void freeOmegastuff(omegastuff *omegas) {

	if (omegas->col) mem_free(omegas->col);
	if (omegas->row) mem_free(omegas->row);
	if (omegas->mean) mem_free(omegas->mean);
	mem_free(omegas);

}//END freeOmegastuff()

void printDecomposeSummary(timeType *tim, probType **prob) {
	int t;

	printf("===============================================================================================================================\n");
	printf("Number of stages                   = %d\n", tim->numStages);
	for ( t = 0; t < tim->numStages; t++ ) {
		printf("-------------------------------------------------------------------------------------------------------------------------------\n");
		printf("Stage %d\n", t);
		printf("Number of decision variables       = %d\t\t", prob[t]->sp->mac);
		printf("(Continuous = %d\tInteger = %d\tBinary = %d)\n", prob[t]->sp->mac - prob[t]->sp->numInt - prob[t]->sp->numBin, prob[t]->sp->numInt, prob[t]->sp->numBin);
		printf("Number of constraints              = %d\n", prob[t]->sp->mar);
		if ( prob[t]->omegas != NULL )
			printf("Number of random variables         = %d\n", prob[t]->omegas->numRV);
		else
			printf("Number of random variables         = 0\n");
	}
	printf("===============================================================================================================================\n");

}//printDecomposeSummary()
