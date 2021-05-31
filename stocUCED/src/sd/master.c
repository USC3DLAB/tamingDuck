/*
 * master.c
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

/* This function is the regularized QP version of master problem. The master problem is solved after the newest cut is added to master problem,
 the incumbent cut is updated if necessary. Here the coefficients on all the cuts are updated, and finally master problem is solved. */
int solveQPMaster(numType *num, sparseVector *dBar, cellType *cell, int IniRow, double lb) {
	double 	d2 = 0.0; /* height at the candidate solution. */
	int 	status, stat1, i;

	if( changeEtaCol(CPXLPptr(cell->master->lp), num->rows, num->cols, cell->k, cell->cuts, lb) ) {
		errMsg("algorithm", "solveMaster", "failed to change the eta column coefficients", 0);
		return 1;
	}

	if ( cell->lbType == NONTRIVIAL ) {
		/* update the right-hand side of cuts to reflect the non-trivial lower bound */
		if ( updateRHS(CPXLPptr(cell->master->lp), cell->cuts, cell->k, cell->lb) ) {
			errMsg("algorithm", "solveQPMaster", "failed to update right-hand side with lower bound information", 0);
			return 1;
		}
	}

#ifdef ALGO_CHECK
	writeProblem(cell->master->lp, "masterCell.lp");
#endif

	/* solve the master problem */
	if ( solveProblem(CPXLPptr(cell->master->lp), cell->master->name, config.MASTERTYPE, &stat1) ) {
		writeProblem(CPXLPptr(cell->master->lp), "error.lp");
		errMsg("algorithm", "solveMaster", "failed to solve the master problem", 0);
		return 1;
	}

	/* increment the number of problems solved during algorithm */
	cell->LPcnt++;

	/* Get the most recent optimal solution to master program */
	status = getPrimal(CPXLPptr(cell->master->lp), cell->candidX, num->cols);
	if ( status ) {
		errMsg("algorithm", "solveMaster", "failed to obtain the primal solution for master", 0);
		return 1;
	}

	/* add the incumbent back to change from \Delta X to X */
	for (i = 1; i <= num->cols; i++)
		d2 += cell->candidX[i] * cell->candidX[i];
	addVectors(cell->candidX, cell->incumbX, NULL, num->cols);

	/* update d_norm_k in soln_type. */
	if (cell->k == 1)
		cell->normDk_1 = d2;
	cell->normDk = d2;

	/* Get the dual solution too */
	status = getDual(CPXLPptr(cell->master->lp), cell->piM, cell->master->mar);
	if ( status ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual solutions to master", 0);
		return 1;
	}
	status = getDualSlacks(CPXLPptr(cell->master->lp), cell->djM, num->cols);
	if ( status ) {
		errMsg("solver", "solveQPMaster", "failed to obtain dual slacks for master", 0);
		return 1;
	}

	/* Find the highest cut at the candidate solution. where cut_height = alpha - beta(xbar + \Delta X) */
	cell->candidEst = vXvSparse(cell->candidX, dBar) + maxCutHeight(cell->cuts, cell->k, cell->candidX, num->cols, lb);

	/* Calculate gamma for next improvement check on incumbent x. */
	cell->gamma = cell->candidEst - cell->incumbEst;

	return 0;
}//END solveQPMaster()

int addCut2Master(cellType *cell, oneCut *cut, int lenX, double lb) {
	intvec 	indices;
	int 	cnt;

	if (!(indices = arr_alloc(lenX + 1, int)))
		errMsg("Allocation", "addcut2Master", "fail to allocate memory to coefficients of beta",0);
	for (cnt = 1; cnt <= lenX; cnt++)
		indices[cnt] = cnt - 1;
	indices[0] = lenX;

	if ( config.MASTERTYPE == PROB_QP )
		cut->alphaIncumb = cut->alpha - vXv(cut->beta, cell->incumbX, NULL, lenX);

	/* check to see if there is room for the candidate cut, else drop a cut */
	if (cell->cuts->cnt == cell->maxCuts) {
		/* make room for the latest cut */
		if( reduceCuts(cell, cell->candidX, cell->piM, lenX, lb) < 0 ) {
			errMsg("algorithm", "addCut2Master", "failed to add reduce cuts to make room for candidate cut", 0);
			return -1;
		}
	}

	/* add the cut to the cell cuts structure as well as on the solver */
	cell->cuts->vals[cell->cuts->cnt] = cut;
	if ( addRow(CPXLPptr(cell->master->lp), lenX + 1, cut->alphaIncumb, GE, 0, indices, cut->beta) ) {
		errMsg("solver", "addcut2Master", "failed to add new row to problem in solver", 0);
		return -1;
	}
	cut->rowNum = cell->master->mar++;

#ifdef CUT_CHECK
	writeProblem(cell->master->lp,"master_wCandidCut.lp");
#endif


	mem_free(indices);
	return cell->cuts->cnt++;
}//END addCuts2Master()

int constructQP(probType *prob, LPptr lp, vectorC incumbX) {
	vectorC rhs;
	intvec indices;
	int cnt;

	if ( !(rhs = (vectorC) arr_alloc(prob->num->rows+1, double)))
		errMsg("allocation", "constructQp", "rhs", 0);
	if ( (!(indices = (intvec) arr_alloc(prob->num->rows, int))) )
		errMsg("allocation", "constructQP", "indices", 0);

	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, incumbX, rhs);

	/* Now we change the right-hand of the master problem. */
	if ( changeRHS(lp, prob->num->rows, indices, rhs + 1) ) {
		errMsg("algorithm", "newMaster", "failed to change the rhs", 0);
		return 1;
	}

	/* change QP bounds */
	if ( changeQPbds(lp, prob->num->cols, prob->sp->bdl, prob->sp->bdu, incumbX) ) {
		errMsg("algorithm", "newMaster", "failed to change the bounds", 0);
		return 1;
	}

	return 0;
}//END constructQP()

/* This function performs the updates on all the coefficients of eta in the master problem constraint matrix.  During every iteration,
 * each of the coefficients on eta are increased, so that the effect of the cut on the objective function is decreased. */
int changeEtaCol(LPptr lp, int numRows, int numCols, int k, cutsType *cuts, double lb) {
	double	etaCoef[1], etaBds[1], coef[1];
	int 	status, c, etaCol[1];
	char	bdsType[1];

	etaCol[0] = numCols;
	bdsType[0] = 'L';

	for (c = 0; c < cuts->cnt; c++){
		/* Currently both incumbent and candidate cuts are treated similarly, and sunk as iterations proceed */
		coef[0] = (double) (k) / (double) cuts->vals[c]->cutObs;         // coeff k/j of eta column

		status = changeCol(lp, numCols, coef, cuts->vals[c]->rowNum, cuts->vals[c]->rowNum+1);
		if ( status ) {
			errMsg("solver", "chgEtaCol", "failed to change eta column in the stage problem", 0);
			return 1;
		}
	}

	if ( cuts->cnt <= 1) {
		if ( cuts->cnt > 0 ) {
			etaCoef[0] = 1.0;
			etaBds[0]  = -INFBOUND;
		}
		else {
			etaCoef[0] = 0.0;
			etaBds[0] = lb;
		}

		status = changeObjx(lp, 1, etaCol, etaCoef);
		if ( status ) {
			errMsg("solver", "changeEtaCol", "failed to change the objective coefficient of eta column in objective function value", 0);
			return 1;
		}

		status = changeBDS(lp, 1, etaCol, bdsType, etaBds);
		if ( status ) {
			errMsg("solver", "changeEtaCol", "failed to change the bound for eta column", 0);
			return 1;
		}
	}

	return 0;
}//END chgEtaCol()

int updateRHS(LPptr lp, cutsType *cuts, int numIter, double lb) {
	int 	cnt;
	vectorC	rhs;
	intvec	indices;

	if (!(rhs = arr_alloc(cuts->cnt, double)))
		errMsg("allocation", "updateRHS", "rhs", 0);
	if (!(indices = arr_alloc(cuts->cnt, int)))
		errMsg("allocation", "updateRHS", "indices", 0);

	for (cnt = 0; cnt < cuts->cnt; cnt++) {
		rhs[cnt] = cuts->vals[cnt]->alphaIncumb + ((double) numIter / (double) cuts->vals[cnt]->cutObs - 1) * lb;
		indices[cnt] = cuts->vals[cnt]->rowNum;
	}

	/* Now we change the right-hand of the master problem. */
	if ( changeRHS(lp, cuts->cnt, indices, rhs) ) {
		errMsg("solver", "changeQPrhs", "failed to change the right-hand side in the solver", 0);
		return 1;
	}

	mem_free(rhs);
	mem_free(indices);

	return 0;
}//END updateRHS

/* Construct the Q diagonal matrix and copy it for quadratic problem. */
int changeQPproximal(LPptr lp, int numCols, double sigma) {
	int    n;
	vectorC qsepvec;

	if (!(qsepvec = arr_alloc(numCols+1, double)))
		errMsg("Allocation", "constructQP", "qsepvec",0);

	/* Construct Q matrix, which is simply a diagonal matrix. */
	for (n = 0; n < numCols; n++)
		qsepvec[n] = 0.5 * sigma;
	qsepvec[n] = 0.0;

	/* Now copy the Q matrix for QP problem. */
	if ( copyQPseparable(lp, qsepvec) ) {
		errMsg("solver", "constructQP", "failed to copy Q matrix", 0);
		return 1;
	}

	mem_free(qsepvec);
	return 0;
}//END constructQP

int changeQPrhs(probType *prob, cellType *cell) {
	int 	cnt, status, offset;
	vectorC 	rhs;
	intvec 	indices;

	if (!(rhs =(vectorC) arr_alloc(cell->master->mar+cell->maxCuts+1, double)))
		errMsg("Allocation", "changeRhs", "rhs",0);
	if (!(indices =(intvec) arr_alloc(cell->master->mar+cell->maxCuts, int)))
		errMsg("Allocation", "changeRhs", "indices",0);

	/* Be careful with the one_norm!! In the CxX() routine, it assumes the 0th element is reserved for the 1_norm, in the returned vector, the T sparse
         vector, and the x vector. */
	for (cnt = 0; cnt < prob->num->rows; cnt++) {
		rhs[cnt + 1] = prob->sp->rhsx[cnt];
		indices[cnt] = cnt;
	}
	offset = prob->num->rows;

	/* b - A * xbar */
	rhs = MSparsexvSub(prob->Dbar, cell->incumbX, rhs);

	/* change the right-hand side of individual cuts for every agent */
	for ( cnt = 0; cnt < cell->cuts->cnt; cnt++ ) {
		rhs[offset + cnt + 1] = cell->cuts->vals[cnt]->alpha - vXv(cell->cuts->vals[cnt]->beta, cell->incumbX, NULL, prob->sp->mac);
		indices[offset + cnt] = cell->cuts->vals[cnt]->rowNum;

		cell->cuts->vals[cnt]->alphaIncumb = rhs[offset + cnt + 1];
	}
	offset += cell->cuts->cnt;

#ifdef RHS_CHECK
	writeProblem(cell->master->lp, "QPchangeRHSBef.lp");
#endif
	/* Now we change the right-hand of the master problem. */
	status = changeRHS(CPXLPptr(cell->master->lp), offset, indices, rhs + 1);
	if ( status ) {
		errMsg("algorithm", "changeRhs", "failed to change the rhs", 0);
		return 1;
	}
#ifdef RHS_CHECK
	writeProblem(cell->master->lp, "QPchangeRHSAft.lp");
#endif

	mem_free(rhs);
	mem_free(indices);
	return 0;
}//END changeQPrh

/* This function changes the (lower) bounds of the variables, while changing from x to d. The lower bounds of d varibles are -xbar
 * (incumbent solution). */
int changeQPbds(LPptr lp, int numCols, vectorC bdl, vectorC bdu, vectorC xk) {
	int 	status = 0, cnt;
	vectorC	lbounds, ubounds;
	intvec	lindices, uindices;
	char 	*llu, *ulu;

	if (!(lbounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "changeBounds", "lbounds",0);
	if (!(lindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "change_bounds", "lindices",0);
	if (!(llu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "llu",0);

	if (!(ubounds = arr_alloc(numCols, double)))
		errMsg("Allocation", "change_bounds", "ubounds",0);
	if (!(uindices = arr_alloc(numCols, int)))
		errMsg("Allocation", "changeBounds", "uindices",0);
	if (!(ulu = arr_alloc(numCols, char)))
		errMsg("Allocation", "changeBounds", "ulu",0);

	/* Change the Upper Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		ubounds[cnt] = bdu[cnt] - xk[cnt + 1];
		uindices[cnt] = cnt;
		ulu[cnt] = 'U';
	}

	status = changeBDS(lp, numCols, uindices, ulu, ubounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the upper bound in the solver", 0);
		return 1;
	}

	/* Change the Lower Bound */
	for (cnt = 0; cnt < numCols; cnt++) {
		lbounds[cnt] = bdl[cnt] - xk[cnt + 1];
		lindices[cnt] = cnt;
		llu[cnt] = 'L';
	}

	status = changeBDS(lp, numCols, lindices, llu, lbounds);
	if (status) {
		errMsg("algorithm", "changeQP", "failed to change the lower bound in the solver", 0);
		return 1;
	}

	mem_free(lbounds); mem_free(lindices); mem_free(llu);
	mem_free(ubounds); mem_free(uindices); mem_free(ulu);

	return 0;
}//END changeQPbds()

/* This subroutine initializes the master problem by copying information from the decomposed prob[0](type: oneProblem) and adding a column for
 * theta for modified benders decomposition. */
oneProblem *newMaster(probType *prob, vectorC xk) {
	oneProblem 	*master;
	int         r, i, j, idx, cnt;
	long        colOffset, rowOffset;
	char        *q;

	if (!(master = (oneProblem *) mem_malloc (sizeof(oneProblem))))
		errMsg("Memory allocation", "new_master", "Faile to allocate memory to mcell->sp", 0);

	/* -+-+-+-+-+-+-+-+-+-+-+-+-+-+- Allocating memory to master -+-+-+-+-+-+-+-+-+-+-+-+-+-+- */
	master->type 	= config.MASTERTYPE;                  	/* type of problem: LP, QP, MIP or MIQP */
	master->objsen 	= prob->sp->objsen;                 	/* sense of the objective: 1 for minimization and -1 for maximization */
	master->mar 	= prob->sp->mar;                       	/* number of rows */
	master->numInt 	= prob->sp->numInt;                 	/* number of integer variables in the problem  */
	master->numnz 	= prob->sp->numnz;                   	/* number of non-zero elements in constraint matrix */
	master->matsz 	= prob->sp->matsz;                   	/* extended matrix size */
	master->marsz 	= prob->sp->marsz;                   	/* extended row size */
	master->rstorsz = prob->sp->rstorsz;               		/* memory size for storing row names */
	master->mac 	= prob->sp->mac+1;           			/* number of columns + etas */
	master->macsz 	= prob->sp->macsz + 1;       			/* extended column size */
	master->cstorsz 	= prob->sp->cstorsz + NAMESIZE;    	/* memory size for storing column names */

	/* Allocate memory to the information whose type is string */
	if (!(master->name = (stringC) arr_alloc(NAMESIZE, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->name",0);
	if (!(master->senx = (stringC) arr_alloc(master->marsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->senx",0);
	if (!(master->ctype = (stringC) arr_alloc(master->macsz,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->ctype",0);
	if (!(master->objname = (stringC) arr_alloc(NAMESIZE,char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objname",0);
	if (!(master->rname = (stringC *) arr_alloc(master->marsz,stringC)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rname",0);
	if (!(master->rstore = (stringC) arr_alloc(master->rstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rstore",0);
	if (!(master->cname = (stringC*) arr_alloc(master->macsz,stringC)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cname",0);
	if (!(master->cstore = (stringC) arr_alloc(master->cstorsz, char)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->cstore",0);

	/* Allocate memory to the information whose type is vector */
	if (!(master->objx = (vectorC) arr_alloc(master->macsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->objx",0);
	if (!(master->rhsx = (vectorC) arr_alloc(master->marsz, double)))
		errMsg("Allocation", "new_master", "Fail to allocate memory to master->rhsx",0);
	if (!(master->matval = (vectorC) arr_alloc(master->matsz, double)))
		errMsg("allocation", "new_master", "master->matval",0);
	if (!(master->bdl = (vectorC) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdl",0);
	if (!(master->bdu = (vectorC) arr_alloc(master->macsz, double)))
		errMsg("allocation", "new_master", "master->bdu",0);

	/* Allocate memory to the information whose type is intvec */
	if (!(master->matbeg = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matbeg",0);
	if (!(master->matcnt = (intvec) arr_alloc(master->macsz, int)))
		errMsg("allocation", "new_master", "master->matcnt",0);
	if (!(master->matind = (intvec) arr_alloc(master->matsz, int)))
		errMsg("allocation", "new_master", "master->matind",0);

	strcpy(master->name, prob->sp->name);           /* Copy problem name */
	strcpy(master->objname, prob->sp->objname);     /* Copy objective name */

	/* Copy problem's column and row names */
	i = 0;
	for (q = prob->sp->cname[0]; q < prob->sp->cname[0] + prob->sp->cstorsz; q++)
		master->cstore[i++] = *q;

	i = 0;
	for (q = prob->sp->rname[0]; q < prob->sp->rname[0] + prob->sp->rstorsz; q++)
		master->rstore[i++] = *q;

	/* Calculate difference in pointers for master/copy row and column names */
	colOffset = master->cstore - prob->sp->cname[0];
	rowOffset = master->rstore - prob->sp->rname[0];

	/* Copy the all column information from the original master problem */
	cnt = 0;
	for (j = 0; j < prob->sp->mac; j++) {
		/* Copy objective function coefficients */
		master->objx[j] = prob->sp->objx[j];
		/* Copy the decision variable type */
		master->ctype[j] = prob->sp->ctype[j];
		/* Copy the upper bound and lower bound */
		master->bdu[j] = prob->sp->bdu[j];
		master->bdl[j] = prob->sp->bdl[j];
		/* Copy column names, offset by length */
		master->cname[j] = prob->sp->cname[j] + colOffset;
		/* Copy the master sparse matrix beginning position of each column */
		master->matbeg[j] = cnt;
		/* Copy the sparse matrix non-zero element count */
		master->matcnt[j] = prob->sp->matcnt[j];
		master->ctype[j] = prob->sp->ctype[j];
		/* Loop through all non-zero elements in this column */
		for (idx = prob->sp->matbeg[j]; idx < prob->sp->matbeg[j] + prob->sp->matcnt[j]; idx++) {
			/* Copy the non-zero coefficient */
			master->matval[cnt] = prob->sp->matval[idx];
			/* Copy the row entry of the non-zero elements */
			master->matind[cnt] = prob->sp->matind[idx];
			cnt++;
		}
	}

	/* Copy all information concerning rows of master */
	for (r = 0; r < prob->sp->mar; r++) {
		/* Copy the right hand side value */
		master->rhsx[r] = prob->sp->rhsx[r];
		/* Copy the constraint sense */
		master->senx[r] = prob->sp->senx[r];
		/* Copy row names, offset by length */
		master->rname[r] = prob->sp->rname[r] + rowOffset;
	}

	/* Initialize information for the extra column in the new master. */
	colOffset = prob->sp->cstorsz;
	strcpy(master->cstore + prob->sp->cstorsz, "eta");
	master->cname[prob->sp->mac] = master->cstore + colOffset;
	master->objx[prob->sp->mac] = 1.0;			// prob->sp->mac is the last column in the original master
	master->ctype[prob->sp->mac] = 'C';
	master->bdu[prob->sp->mac] = INFBOUND;
	master->bdl[prob->sp->mac] = 0.0;
	master->matbeg[prob->sp->mac] = prob->sp->numnz;	// Beginning point in matval/matind in eta columns. every eta column begins at the same address
	master->matcnt[prob->sp->mac] = 0;               // Only optimality cuts has eta

	/* Load the copy into CPLEX */
	master->lp = setupProblem(master->name, master->type, master->mac, master->mar, master->objsen, master->objx, master->rhsx, master->senx, master->matbeg, master->matcnt,master->matind, master->matval, master->bdl, master->bdu, NULL, master->cname, master->rname, master->ctype);
	if ( master->lp == NULL ) {
		errMsg("Problem Setup", "new_master", "failed to setup master problem in the solver",0);
		return NULL;
	}

#if 0
	status = writeProblem(master->lp, "newMaster.lp");
	if ( status ) {
		errMsg("write problem", "new_master", "failed to write master problem to file",0);
		return NULL;
	}
#endif

	return master;

}//END newMaster
