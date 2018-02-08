/*
 * input.c
 *
 *  Created on: Sep 29, 2015
 *      Author: gjharsha
 */

#include <math.h>
#include "smps.h"
#include "solver.h"
#include <stdio.h>
#include <string.h>
#include "utils.h"

int readFiles(stringC inputDir, stringC probName, oneProblem **orig, timeType **tim, stocType **stoc) {

	/* read problem core file */
	(*orig) = readCore(inputDir, probName);
	if ( (*orig) == NULL ) {
		errMsg("read", "readFiles", "failed to read problem core file", 0);
		return 1;
	}

	/* read problem time file */
	(*tim) = readTime(inputDir, probName, (*orig));
	if ( (*tim) == NULL ) {
		errMsg("read", "readFiles", "failed to read problem time file", 0);
		return 1;
	}

	(*stoc) = readStoc(inputDir, probName, (*orig), (*tim));
	if ( (*stoc) == NULL ) {
		errMsg("read", "readFiles", "failed to read problem stoc file", 0);
		return 1;
	}

	printf("Successfully read all the files for problem: %s.\n\n", probName);

#ifdef INPUT_CHECK
	writeStocType((*stoc));
#endif

	return 0;
}//END readFiles()

oneProblem *readCore(stringC inputDir, stringC probName) {
	LPptr 			lp = NULL;
	char 			probpath[BLOCKSIZE], line[BLOCKSIZE], field1[NAMESIZE], field2[NAMESIZE];
	oneProblem      *orig;
	int				c, nzcnt;
	FILE			*fptr;

	/* Locate the problem core file */
	sprintf(probpath, "%s%s/%s.cor", inputDir, probName, probName);

	fptr = fopen(probpath, "r");
	if ( fptr == NULL ) {
		sprintf(probpath, "%s%s/%s.mps", inputDir, probName, probName);
		fptr = fopen(probpath, "r");
		if ( fptr == NULL ) {
			errMsg("read", "readCore", "failed to open problem core file", 0);
			return NULL;
		}
	}

	/* NAME section: read problem name and compare with that read earlier */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", field1, field2);
	else {
		errMsg("read", "readCore", "failed to read the problem name", 0);
		return NULL;
	}
	if ( !(strncmp(field1, "NAME", 4)) )
		if( strcmp(probName, field2) ) {
			errMsg("read", "readCore", "problem name does not match, in NAME section", 0);
			return NULL;
		}
	fclose (fptr);
	/* Create LP pointer */
	if ((createProblem(probName, &lp))) {
		errMsg("solver", "readCore", "failed to create problem in solver.\n", 0);
		return NULL;
	}

	if ((readProblem(probpath, lp))) {
		errMsg("solver", "readCore", "failed to read and copy the problem data", 0);
		return NULL;
	}

	/* Allocate memory to the elements of problem and assign default values*/
	orig = (oneProblem *) mem_malloc(sizeof(oneProblem));
	orig->lp = lp;
	orig->name = (stringC) mem_calloc(NAMESIZE, sizeof(char));

	/* obtain type of problem read */
	orig->type = getProbType(lp);

	/* obtain the problem elements */
	/* (1) problem name */
	strcpy(orig->name, probName);

	/* (2) objective sense */
	orig->objsen = getObjSen(lp);
	if (!(orig->objsen))
		errMsg("solver", "readCore", "failed to obtain the objective sense", 1);

	/* (3) number of rows */
	orig->mar = getNumRows(lp);
	if (!(orig->mar))
		errMsg("solver", "readCore", "failed to obtain the number of rows in the problem", 1);

	/* (4) number of columns */
	orig->mac = getNumCols(lp);
	if (!(orig->mac))
		errMsg("solver", "readCore", "failed to obtain the number of columns in the problem", 1);

	/* (5) number of non-zeros */
	nzcnt = getNumnz(lp);

	/* continue allocating memory to the elements of problem and assign default values*/
	orig->objx = (vectorC) mem_calloc(orig->mac, sizeof(double));
	orig->rhsx = (vectorC) mem_calloc(orig->mar, sizeof(double));
	orig->senx = (stringC) mem_malloc(orig->mar*sizeof(char));
	orig->matbeg = (intvec) mem_malloc(orig->mac*sizeof(int));
	orig->matcnt = (intvec) mem_malloc(orig->mac*sizeof(int));
	orig->matind = (intvec) mem_malloc(nzcnt*sizeof(int));
	orig->matval = (vectorC) mem_malloc(nzcnt*sizeof(double));
	orig->bdl = (vectorC) mem_malloc(orig->mac*sizeof(double));
	orig->bdu = (vectorC) mem_malloc(orig->mac*sizeof(double));
	orig->ctype = (stringC) mem_malloc(orig->mac*sizeof(char));

	/* (6) objective function coefficients */
	if ( (getObjx(lp, 0, orig->mac, orig->objx)) )
		errMsg("solver", "readCore", "failed to obtain the objective coefficients", 1);

	/* (7) constraint right hand side */
	if ( (getRhsx(lp, 0, orig->mar, orig->rhsx)) )
		errMsg("solver", "readCore", "failed to obtain problem constraint right hand sides", 1);

	/* (8) constraint sense */
	if ( (getSense(lp, 0, orig->mar, orig->senx)))
		errMsg("solver", "readCore", "failed to obtain problem constraint sense", 1);

	/* (10) constraint matrix coefficients */
	if ( (getCols(lp, 0, orig->mac, orig->matbeg, orig->matind, orig->matval, nzcnt)) )
		errMsg("solver", "readCore", "failed to obtain the constraint matrix coefficients", 1);
	for ( c = 0; c < orig->mac - 1; c++ )
		orig->matcnt[c] = orig->matbeg[c+1] - orig->matbeg[c];
	orig->matcnt[c] = nzcnt - orig->matbeg[c];

	/* (11) problem variable bounds */
	if( getLb(lp, 0, orig->mac, orig->bdl) )
		errMsg("solver", "readCore", "failed to obtain the problem lower bounds", 1);
	if( getUb(lp, 0, orig->mac, orig->bdu) )
		errMsg("solver", "readCore", "failed to obtain the problem upper bounds", 1);

	if ( orig->type == PROB_MILP || orig->type == PROB_MIQP ) {
		/* (12) get problem constraint type */
		if ( getCtype(lp, 0, orig->mac, orig->ctype) )
			errMsg("solver", "readCore", "failed to obtain the variable types", 1);
		orig->numInt = getNumInt(lp);
		orig->numBin = getNumBinary(lp);
	}
	else {
		for ( c = 0; c < orig->mac; c++ )
			orig->ctype[c] = 'C';
		orig->numInt = 0;
	}

	/* Allocate memory to hold the names of problem elements */
	orig->objname = (stringC) mem_calloc(NAMESIZE, sizeof(char));
	orig->cstorsz = -getCstoreSize(lp, 0, orig->mac);
	if ( orig->cstorsz <= 0 )
		errMsg("solver", "readCore", "Could not determine amount of space for column names", 1);
	orig->cname = (stringC *) mem_malloc(orig->mac*sizeof(char *));
	orig->cstore = (stringC) mem_malloc(orig->cstorsz);

	orig->rstorsz = -getRstoreSize(lp, 0, orig->mar);
	if ( orig->rstorsz < 0 )
		errMsg("solver", "readCore", "Could not determine amount of space for row names", 1);
	orig->rname = (stringC *) mem_malloc(orig->mar*sizeof(char *));
	orig->rstore = (stringC) mem_malloc(orig->rstorsz);

	/* (12) objective name */
	if ( (getObjName(lp, orig->objname)) )
		errMsg("solver", "readCore", "failed to obtain objective name", 1);

	/* (13) problem row name */
	if ( (getRowName(lp, 0, orig->mar, orig->rname, orig->rstore, orig->rstorsz)) )
		errMsg( "solver", "readCore", "failed to obtain row names", 1);

	/* (14) problem column name */
	if ( (getColName(lp, 0, orig->mac, orig->cname, orig->cstore, orig->cstorsz)) )
		errMsg("solver", "readCore", "failed to obtain column names", 1);

	orig->matsz = nzcnt;
	orig->macsz = nzcnt;
	orig->marsz = nzcnt;
	orig->numnz = nzcnt;

	return orig;

}//END readCore()

timeType *readTime(stringC inputDir, stringC probName, oneProblem *orig) {
	timeType	*tim;
	char		probpath[2*BLOCKSIZE], line[BLOCKSIZE], field1[NAMESIZE], field2[NAMESIZE];
	int			defaultStages = 20, n, m;
	FILE		*fptr;

	/* Locate the problem core file */
	sprintf(probpath, "%s%s/%s.tim", inputDir, probName, probName);

	/* open the time file */
	fptr = fopen(probpath,"r");
	if (fptr == NULL) {
		errMsg("read","readTime","failed to read time file/missing file", 0);
		return NULL;
	}

	/* allocate memory and initialize */
	if (!(tim = (timeType *) mem_malloc(sizeof(timeType))))
		errMsg("allocation", "readTime", "timeType",0);
	if(!(tim->stgNames = (stringC *) mem_malloc(defaultStages*sizeof(stringC))))
		errMsg("allocation", "readTime", "stgNames in timeType", 0);
	tim->numStages = 0; n = 0;
	tim->numCols = 0; tim->numRows = 0;

	/* TIME section: collect the problem name */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", field1, field2);
	else {
		errMsg("read", "readTime", "failed to read the problem name", 0);
		return NULL;
	}
	if ( !(strncmp(field1, "TIME", 4)) )
		if( strcmp(probName, field2) ) {
			errMsg("read", "readTime", "problem name does not match, TIME section", 0);
			return NULL;
		}

	/* PERIODS section: collect the time file type */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", field1, field2);
	else {
		errMsg("read", "readTime", "failed to read the time file type", 0);
		return NULL;
	}

	if ( !(strncmp(field1, "PERIODS", 7)) ) {
		if ( !(strncmp(field2, "EXPLICIT", 8)) )
			tim->type = 1;
		else
			tim->type = 0;
	}
	else {
		errMsg("read", "readTime", "unknown header record time file, PERIODS section", 0);
		return NULL;
	}

	if ( tim->type == 0 ){
		if( !(tim->row = (intvec) arr_alloc(defaultStages, int)) )
			errMsg("allocation", "readTime", "rowNames in timeType", 0);
		if( !(tim->col = (intvec) arr_alloc(defaultStages, int)) )
			errMsg("allocation", "readTime", "colNames in timeType", 0);
		while ( fgets(line, sizeof line, fptr )!= NULL ) {
			if (strncmp(line,"ENDATA",6)) {
				if ( !(tim->stgNames[n] =  (stringC) mem_malloc(NAMESIZE*sizeof(char))))
					errMsg("allocation", "readTime", "individual stage names", 0);
				sscanf(line, "%s %s %s", field1, field2, tim->stgNames[n]);
				/* find the column and row coordinates in original problem */
				m = 0;
				while ( m < orig->mac ){
					if ( !(strcmp(field1, orig->cname[m])) )
						break;
					m++;
				}
				if ( m == orig->mac ) {
					errMsg("read", "readTime", "unknown column name in the time file", 0);
					return NULL;
				}
				tim->col[n] = m;
				if ( !(strcmp(field2, orig->objname)) )
					m = 0;
				else {
					m = 0;
					while (m < orig->mar ) {
						if ( !(strcmp(field2, orig->rname[m])) )
							break;
						m++;
					}
				}
				if ( m == orig->mar ) {
					errMsg("read", "readTime", "unknown row name in the time file", 0);
					return NULL;
				}
				tim->row[n] = m;

				/* increment stage counter and proceed */
				tim->numStages++; n++;
				if ( n == defaultStages ) {
					errMsg("allocation", "readTime",
							"ran out of memory to store periods in implicit declaration; increase default value", 0);
					return NULL;
				}
			}
		}
		tim->numCols = tim->numRows = tim->numStages;
		tim->rowStg = NULL; tim->colStg = NULL;
	}
	else {
		errMsg("read", "readTime", "explicit time file description is currently not supported", 0);
		/* TODO (HG): complete the code to read explicit time file declaration */
		return NULL;
	}

	/* reallocate memory elements of time structure */
	tim->stgNames = (stringC *) mem_realloc(tim->stgNames, tim->numStages*sizeof(stringC));
	if (tim->type == 0) {
		tim->col = (intvec) mem_realloc(tim->col, tim->numStages*sizeof(int));
		tim->row = (intvec) mem_realloc(tim->row, tim->numStages*sizeof(int));
	}
	else {
		tim->col = (intvec) mem_realloc(tim->col, tim->numStages*sizeof(int));
		tim->row = (intvec) mem_realloc(tim->row, tim->numStages*sizeof(int));
		tim->colStg = (intvec) mem_realloc(tim->colStg, tim->numCols*sizeof(int));
		tim->rowStg = (intvec) mem_realloc(tim->rowStg, tim->numRows*sizeof(int));
	}
	fclose(fptr);
	return tim;
}//END readTime()

stocType *readStoc(stringC inputDir, stringC probName, oneProblem *orig, timeType *tim) {
	stocType *stoc;
	char	probpath[2*BLOCKSIZE], line[BLOCKSIZE], **fields, fieldType;
	FILE	*fptr;
	int		maxOmegas = 1000, maxVals = 4000, n, numFields, maxFields = 10;

	/* Locate the problem sto file */
	sprintf(probpath, "%s%s/%s.sto", inputDir, probName, probName);

	/* open the time file */
	fptr = fopen(probpath,"r");
	if (fptr == NULL) {
		errMsg("read","readStoch","failed to read stoch file file/missing file", 0);
		return NULL;
	}

	/* allocate memory to field locations */
	if (!(fields = (stringC *) arr_alloc(maxFields, stringC)) )
		errMsg("allocation", "readStoc", "field locations", 0);
	for (n = 0; n < maxFields; n++ )
		if ( !(fields[n] = (stringC) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "readStoc", "individual field location", 0);

	/* allocate memory to stocType and initialize elements */
	if ( !(stoc = (stocType *) mem_malloc(sizeof(stocType))) )
		errMsg("allocation", "readStoc", "stoc", 0);
	if ( !(stoc->col = (intvec) arr_alloc(maxOmegas, int)) )
		errMsg("allocation", "readStoc", "stoc->col", 0);
	if ( !(stoc->row = (intvec) arr_alloc(maxOmegas, int)) )
		errMsg("allocation", "readStoc", "stoc->row", 0);
	if ( !(stoc->mean = (vectorC) arr_alloc(maxOmegas, double)) )
		errMsg("allocation", "readStoc", "stoc->mean", 0);
	if ( !(stoc->numVals = (intvec) arr_alloc(maxOmegas, int)) )
		errMsg("allocation", "readStoc", "stoc->numVals", 0);
	if ( !(stoc->vals = (vectorC *) arr_alloc(maxOmegas, vectorC)) )
		errMsg("allocation", "readStoc", "stoc->vals", 0);
	if ( !(stoc->probs = (vectorC *) arr_alloc(maxOmegas, vectorC)) )
		errMsg("allocation", "readStoc", "stoc->vals", 0);
	if ( !(stoc->type = (stringC) arr_alloc(NAMESIZE, char)) )
		errMsg("allocation", "readStoc", "stoc->type", 0);
	if ( !(stoc->groupBeg = (intvec) arr_alloc(maxOmegas, int)) )
		errMsg("allocation", "readStoc", "stoc->groupBeg", 0);
	if ( !(stoc->numPerGroup = (intvec) arr_alloc(maxOmegas, int)) )
		errMsg("allocation", "readStoc", "stoc->numPerGroup", 0);
	stoc->numCipher = 0;
	stoc->numOmega = 0;
	stoc->numGroups = 0;
	stoc->sim = CFALSE;
	stoc->arma = NULL;

	/* STOCH section: read problem name and compare with that read earlier */
	if ( fgets(line, sizeof line, fptr) != NULL )
		sscanf(line, "%s %s", fields[0], fields[1]);
	else {
		errMsg("read", "readStoc", "failed to read the problem name", 0);
		return NULL;
	}
	if ( !(strncmp(fields[0], "STOCH", 5)) )
		if( strcmp(probName, fields[1]) ) {
			errMsg("read", "readStoc", "Problem name do not match, in STOCH section", 0);
			return NULL;
		}

	while ( !(getLine(&fptr, fields, &fieldType, &numFields)) ) {
		if ( !(strcmp(fields[0], "INDEP")) ) {
			if ( readIndep(fptr, fields, orig, maxOmegas, maxVals, stoc) ) {
				errMsg("read", "readStoc", "failed to read stoch file with independent data", 0);
				return NULL;
			}
		}
		else if ( !(strcmp(fields[0], "BLOCKS")) ) {
			if ( readBlocks(fptr, fields, orig, maxOmegas, maxVals, stoc) ) {
				errMsg("read", "readStoc", "failed to read stoch file with blocks", 0);
				return NULL;
			}
		}
		else if ( !(strcmp(fields[0], "SCENARIOS")) ) {
			if ( readScenarios(fptr, fields, orig, tim, maxOmegas, maxVals, stoc) ) {
				errMsg("read", "readStoc", "failed to read stoch file with scenarios", 0);
				return NULL;
			}
		}
		else if ( !(strcmp(fields[0], "ENDATA")) )
			break;
		else
			continue;
	}

	/* free allocated memory */
	for (n = 0; n < maxFields; n++ )
		if (fields[n]) mem_free(fields[n]);
	mem_free(fields);
	fclose(fptr);

	return stoc;
}//END readStoc()

int readIndep(FILE *fptr, stringC *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc) {
	stringC 	*rvRows, *rvCols;
	char	strType;
	int		n, numFields;

	/* allocate memory to hold the names of random variable */
	if ( !(rvRows = (stringC *) arr_alloc(maxOmegas, stringC)) )
		errMsg("allocation", "readIndep", "rvNames", 0);
	if ( !(rvCols = (stringC *) arr_alloc(maxOmegas, stringC)) )
		errMsg("allocation", "readIndep", "rvNames", 0);

	if ( !(strcmp(fields[1], "DISCRETE")) ) {
		/* store the type of stochastic process encountered */
		sprintf(stoc->type, "INDEP_DISCRETE");

		while (CTRUE) {
			getLine(&fptr, fields, &strType, &numFields);
			if (strType != 'f')
				break;										//Encountered ENDATA
			n = stoc->numOmega - 1;
			if ( n > maxOmegas ) {
				errMsg("allocation", "readIndep", "reached maxOmega limit for INDEP format", 0);
				return 1;
			}
			while (n >= 0 ) {
				if ( !(strcmp(fields[0], rvCols[n])) && !(strcmp(fields[1], rvRows[n])) )
					break;
				n--;
			}
			if ( n == -1 ) {
				/* new random variable encountered */
				if ( !(rvRows[stoc->numOmega] = (stringC) arr_alloc(NAMESIZE, char)) )
					errMsg("allocation", "readIndep", "rvNames[n]", 0);
				if ( !(rvCols[stoc->numOmega] = (stringC) arr_alloc(NAMESIZE, char)) )
					errMsg("allocation", "readIndep", "rvNames[n]", 0);
				if ( !(stoc->vals[stoc->numOmega] = (vectorC) arr_alloc(maxVals, double)) )
					errMsg("allocation", "readIndep","omega.vals[n]", 0);
				if ( !(stoc->probs[stoc->numOmega] = (vectorC) arr_alloc(maxVals, double)) )
					errMsg("allocation", "readIndep", "omega.probs[n]", 0);

				strcpy(rvCols[stoc->numOmega], fields[0]);
				strcpy(rvRows[stoc->numOmega], fields[1]);
				stoc->numVals[stoc->numOmega++] = 0;
				/* identify row and column coordinates in the problem */
				if ( !(strcmp(fields[0], "RHS")) )
					n = -1;
				else {
					n = 0;
					while ( n < orig->mac ){
						if ( !(strcmp(rvCols[stoc->numOmega-1], orig->cname[n])) )
							break;
						n++;
					}
				}
				if ( n == orig->mac ) {
					errMsg("read", "readIndep", "unknown column name in the stoch file", 0);
					return 1;
				}
				stoc->col[stoc->numOmega-1] = n;
				if ( !(strcmp(fields[1], orig->objname)) )
					n = -1;
				else {
					n = 0;
					while (n < orig->mar ) {
						if ( !(strcmp(rvRows[stoc->numOmega-1], orig->rname[n])) )
							break;
						n++;
					}
				}
				if ( n == orig->mar ) {
					errMsg("read", "readIndep", "unknown row name in the stoch file", 0);
					return 1;
				}
				stoc->row[stoc->numOmega-1] = n;
			}
			if ( numFields == 4) {
				stoc->vals[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[2]);
				stoc->probs[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[3]);
				stoc->mean[stoc->numOmega-1] += str2float(fields[2])*str2float(fields[3]);
				stoc->numVals[stoc->numOmega-1]++;
			}
			else if ( numFields == 5 ) {
				stoc->vals[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[2]);
				stoc->probs[stoc->numOmega-1][stoc->numVals[stoc->numOmega-1]] = str2float(fields[4]);
				stoc->mean[stoc->numOmega-1] += str2float(fields[2])*str2float(fields[4]);
				stoc->numVals[stoc->numOmega-1]++;
			}
			else {
				errMsg("read", "readIndep", "missing field in stoch file", 0);
				return 1;
			}
		}
	}
	else if ( strstr(fields[1], "NORMAL") != NULL ) {
		/* store the type of stochastic process encountered */
		if ( strcmp(fields[1], "NORMAL") )
			stoc->sim = CTRUE;
		sprintf(stoc->type, "INDEP_%s",fields[1]);

		if ( !(stoc->vals[0] = (vectorC) arr_alloc(maxOmegas, double)) )
			errMsg("allocation", "readIndep","omega.vals[n]", 0);
		mem_free(stoc->probs); stoc->probs = NULL;

		while (CTRUE) {
			getLine(&fptr, fields, &strType, &numFields);
			if (strType != 'f')
				break;
			n = stoc->numOmega - 1;
			if ( n > maxOmegas ) {
				errMsg("allocation", "readIndep", "reached maxOmega limit for INDEP format", 0);
				return 1;
			}
			while (n >= 0 ) {
				if ( !(strcmp(fields[0], rvCols[n])) && !(strcmp(fields[1], rvRows[n])) )
					break;
				n--;
			}
			if ( n == -1 ) {
				/* new random variable encountered */
				if ( stoc->numOmega == maxOmegas ) {
					errMsg("read", "readIndep", "ran out of memory to store row and column names", 0);
					return 1;
				}
				if ( !(rvRows[stoc->numOmega] = (stringC) arr_alloc(NAMESIZE, char)) )
					errMsg("allocation", "readIndep", "rvNames[n]", 0);
				if ( !(rvCols[stoc->numOmega] = (stringC) arr_alloc(NAMESIZE, char)) )
					errMsg("allocation", "readIndep", "rvNames[n]", 0);

				strcpy(rvCols[stoc->numOmega], fields[0]);
				strcpy(rvRows[stoc->numOmega], fields[1]);
				stoc->numVals[stoc->numOmega++] = 0;
				/* identify row and column coordinates in the problem */
				if ( !(strcmp(fields[0], "RHS")) )
					n = -1;
				else {
					n = 0;
					while ( n < orig->mac ){
						if ( !(strcmp(rvCols[stoc->numOmega-1], orig->cname[n])) )
							break;
						n++;
					}
				}
				if ( n == orig->mac ) {
					errMsg("read", "readIndep", "unknown column name in the stoch file", 0);
					return 1;
				}
				stoc->col[stoc->numOmega-1] = n;
				if ( !(strcmp(fields[1], orig->objname)) )
					n = -1;
				else {
					n = 0;
					while (n < orig->mar ) {
						if ( !(strcmp(rvRows[stoc->numOmega-1], orig->rname[n])) )
							break;
						n++;
					}
				}
				if ( n == orig->mar ) {
					errMsg("read", "readIndep", "unknown row name in the stoch file", 0);
					return 1;
				}
				stoc->row[stoc->numOmega-1] = n;
			}
			if ( numFields == 4) {
				/* note, standard deviation is held in the first vals field */
				stoc->mean[stoc->numOmega-1] 	= str2float(fields[2]);
				stoc->vals[0][stoc->numOmega-1] = sqrt(str2float(fields[3]));
			}
			else if ( numFields == 5 ) {
				stoc->mean[stoc->numOmega-1] 	= str2float(fields[2]);
				stoc->vals[0][stoc->numOmega-1] = sqrt(str2float(fields[4]));
			}
			else {
				errMsg("read", "readIndep", "missing field in stoch file", 0);
				return 1;
			}
		}
	}
	else if ( !(strcmp(fields[1], "EXPONENTIAL")) ) {
		errMsg("read", "readIndeps", "no support for exponential distribution type in INDEP section", 1);
		return 1;
	}
	else if ( !(strcmp(fields[1], "UNIFORM")) ) {
		errMsg("read", "readIndeps", "no support for uniform distribution type in INDEP section", 1);
		return 1;
	}
	else if ( !(strcmp(fields[1], "GAMMA")) ) {
		errMsg("read", "readIndeps", "no support for gamma distribution type in INDEP section", 1);
		return 1;
	}
	else if ( !(strcmp(fields[1], "GEOMETRIC")) ) {
		errMsg("read", "readIndeps", "no support for geometric distribution type in INDEP section", 1);
		return 1;
	}
	else {
		errMsg("read", "readIndeps", "unknown distribution type in INDEP section", 1);
		return 1;
	}

	for ( n = 0; n < stoc->numOmega; n++ ) {
		mem_free(rvCols[n]); mem_free(rvRows[n]);
	}
	mem_free(rvCols); mem_free(rvRows);

	/* increase the number of stochastic variables groups */
	stoc->groupBeg[stoc->numGroups] = 0;
	stoc->numPerGroup[stoc->numGroups++] = stoc->numOmega;

	return 0;
}//END readIndep()

int readBlocks(FILE *fptr, stringC *fields, oneProblem *orig, int maxOmegas, int maxVals, stocType *stoc) {
	int status;

	if ( !(strcmp(fields[1], "DISCRETE")) ) {
		/* store the type of stochastic process encountered */
		sprintf(stoc->type, "BLOCKS_DISCRETE");
		status = readBlk(fptr, fields, orig, maxOmegas, maxVals, CTRUE, stoc);
		if ( status ) {
			errMsg("read", "readBlocks", "failed to read independent blocks structure", 0);
			return 1;
		}
	}
	else if ( !(strcmp(fields[1], "MVNORMAL")) ) {
		errMsg("read", "readBlocks", "no support for multivariate normal distribution type in BLOCKS section", 0);
		return 1;
	}
	else if ( !(strcmp(fields[1], "LINTRAN")) ) {
		errMsg("read", "readBlocks", "no support for linear translation type in BLOCKS section", 1);
		return 1;
	}
	else if ( !(strcmp(fields[1], "ARMA")) ) {
		status = readARMA(fptr, fields, orig, stoc, maxOmegas);
		if ( status ) {
			errMsg("read", "readBlocks", "failed to read ARMA blocks structure", 0);
			return 1;
		}
	}
	else {
		errMsg("read", "readBlocks", "unknown distribution type in BLOCKS section", 1);
		return 1;
	}

	return 0;
}//END readBlocks()

int readBlk(FILE *fptr, stringC *fields, oneProblem *orig, int maxOmegas, int maxVals, BOOL origRV, stocType *stoc) {
	stringC 	*rvRows, *rvCols;
	char 	strType, currBlock[NAMESIZE] = "\0";
	int		numFields, numRV=0, n;
	BOOL	newBlk;

	/* allocate memory to hold the names of random variable */
	if ( !(rvRows = (stringC *) arr_alloc(maxOmegas, stringC)) )
		errMsg("allocation", "readIndep", "rvNames", 0);
	if ( !(rvCols = (stringC *) arr_alloc(maxOmegas, stringC)) )
		errMsg("allocation", "readIndep", "rvNames", 0);

	for ( n = 0; n < maxOmegas; n++) {
		if ( !(rvCols[n] = (stringC) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "readBlk", "rvCols", 0);
		if ( !(rvRows[n] = (stringC) arr_alloc(NAMESIZE, char)) )
			errMsg("allocation", "readBlk", "rvRows", 0);
	}

	while (CTRUE) {
		getLine(&fptr, fields, &strType, &numFields);
		if (strType != 'f')
			break;
		if ( !(strcmp(fields[0], "BL")) ) {
			/* new realization of the block */
			if ( strcmp(currBlock, fields[1]) ) {
				/* first encounter with the block, prepare to record names of random variables */
				newBlk = CTRUE;
				strcpy(currBlock, fields[1]);
				stoc->groupBeg[stoc->numGroups] = stoc->numOmega;
				stoc->numPerGroup[stoc->numGroups] = numRV = 0;
				if ( !(stoc->probs[stoc->numGroups] = (vectorC) arr_alloc(maxVals, double)) )
					errMsg("allocation", "readBlk", "stoc->prob[n]", 0);
				stoc->probs[stoc->numGroups][stoc->numVals[stoc->numGroups]++] = str2float(fields[3]);
				stoc->numGroups++;
			}
			else {
				newBlk = CFALSE;
				if ( stoc->numVals[stoc->numGroups-1] == maxVals )
					errMsg("allocation", "readBlock", "exceeded memory limit on maxVals", 1);
				stoc->probs[stoc->numGroups-1][stoc->numVals[stoc->numGroups-1]++] = str2float(fields[3]);
			}
		}
		else {
			/* read block elements */
			if (newBlk) {
				/* record names of random variables on their first realization */
				strcpy(rvCols[numRV], fields[0]);
				if ( origRV ) {
					/* if the random variables are problem random variables then identify their coordinates in original problem */
					strcpy(rvRows[numRV], fields[1]);
				}
				/* column coordinates */
				if ( !(strcmp(fields[0], "RHS")) )
					n = -1;
				else {
					n = 0;
					while ( n < orig->mac ) {
						if ( !(strcmp(fields[0], orig->cname[n])) )
							break;
						n++;
					}
				}
				stoc->col[stoc->numOmega] = n;
				/* row coordinates */
				if ( !(strcmp(fields[1], orig->objname)) )
					n = -1;
				else {
					n = 0;
					while ( n < orig->mar ) {
						if ( !(strcmp(fields[1], orig->rname[n])) )
							break;
						n++;
					}
				}
				stoc->row[stoc->numOmega] = n;

				/* make sure there is memory space available for new realization and store it */
				if (stoc->numOmega == maxOmegas )
					errMsg("allocation", "readBlock", "reached max limit maxOmegas", 1);

				if ( !(stoc->vals[stoc->numOmega] = (vectorC) arr_alloc(maxVals, double)) )
					errMsg("allocation", "readBlock","omega.vals[n]", 0);

				if (origRV == 1) {
					stoc->vals[stoc->numOmega][stoc->numVals[stoc->numGroups-1]-1] = str2float(fields[2]);
					stoc->mean[stoc->numOmega] += stoc->probs[stoc->numGroups-1][stoc->numVals[stoc->numGroups-1]-1]*str2float(fields[2]);
					stoc->numOmega++;
				}
				else {
					errMsg("read", "readBlk", "reading auxiliary variables not supported", 0);
					return 1;
				}
				/* increment the number of random variables in the group */
				stoc->numPerGroup[stoc->numGroups-1]++;
				numRV++;
			}
			else {
				/* locate the random variable in the list and record realization */
				n = 0;
				while (n < numRV) {
					if ( origRV == 	CFALSE && !(strcmp(rvCols[n], fields[0])) )
						break;
					else if ( origRV == CTRUE && !(strcmp(rvCols[n], fields[0])) && !(strcmp(rvRows[n], fields[1])) )
						break;
					n++;
				}
				if ( n == numRV )
					errMsg("read", "readBlock", "unknown block random variable name", 1);

				n += stoc->groupBeg[stoc->numGroups-1];
				if ( origRV == 1 ) {
					/* the third field has values */
					stoc->vals[n][stoc->numVals[stoc->numGroups-1]-1] = str2float(fields[2]);
					stoc->mean[n] += stoc->probs[stoc->numGroups-1][stoc->numVals[stoc->numGroups-1]-1]*str2float(fields[2]);
				}
				else {
					/* the second field has values */
					errMsg("read", "readBlk", "reading auxiliary variables not supported", 0);
					return 1;
				}
			}
		}
	}

	for ( n = 0; n < maxOmegas; n++) {
		if ( rvRows[n] ) mem_free(rvRows[n]);
		if ( rvCols[n] ) mem_free(rvCols[n]);
	}
	mem_free(rvRows); mem_free(rvCols);

	return 0;
}//END readBlk()

int readARMA(FILE *fptr, stringC *fields, oneProblem *orig, stocType *stoc, int maxOmegas) {
	armaType *arma;
	char 	strType, currBlock[NAMESIZE] = "\0";
	int		numFields, numRV = 0, maxP = 10, maxQ = 10, maxEps = 10, maxT = 365, armaR = 0, armaE = 0, armaS = 0, j, n;
	BOOL	newBlk = CFALSE;

	/* set up elements of stocType. Free-up elements not used, and trim down the number of groups to two (one for problem stochastic elements and
	 * one for the residuals) */
	strcpy(stoc->type, "ARMA");
	stoc->sim = CFALSE;
	mem_free(stoc->numVals); stoc->numVals = NULL;
	mem_free(stoc->vals); stoc->vals = NULL;
	mem_free(stoc->probs); stoc->probs = NULL;
	stoc->numPerGroup = (intvec) mem_realloc(stoc->numPerGroup, 2*sizeof(int));
	stoc->groupBeg = (intvec) mem_realloc(stoc->groupBeg, 2*sizeof(int));

	/* allocate memory to hold information about the ARMA process */
	if ( !(arma = (armaType *) mem_malloc(sizeof(armaType))) )
		errMsg("allocation", "readARMA", "armaType", 0);
	arma->obs = arma->eps = arma->eta = arma->sigma = NULL;
	arma->AR = arma->MA = NULL; arma->meanEps = arma->varEps = NULL;
	arma->p = arma->q = 0; arma->T = arma->N = 0;

	while (CTRUE) {
		getLine(&fptr, fields, &strType, &numFields);
		if (strType != 'f')
			break;
		if ( !(strcmp(fields[0], "BL")) ) {
			/* new block encountered */
			strcpy(currBlock, fields[1]);
			newBlk = CTRUE;
		}
		else {
			if ( !(strcmp(currBlock, "VAR")) ) {
				/* stochastic elements in the problem */
				/* column coordinates */
				if ( !(strcmp(fields[0], "RHS")) )
					j = -1;
				else {
					j = 0;
					while ( j < orig->mac ) {
						if ( !(strcmp(fields[0], orig->cname[j])) )
							break;
						j++;
					}
				}
				stoc->col[stoc->numOmega] = j;
				if ( !(strcmp(fields[1], orig->objname)) )
					j = -1;
				else {
					j = 0;
					while ( j < orig->mar ) {
						if ( !(strcmp(fields[1], orig->rname[j])) )
							break;
						j++;
					}
				}
				stoc->row[stoc->numOmega] = j;
				stoc->numOmega++; stoc->numPerGroup[0]++;

				/* make sure there is memory space available for new realization and store it */
				if (stoc->numOmega == maxOmegas )
					errMsg("allocation", "readBlock", "reached max limit maxOmegas", 1);

			}
			else if ( !(strcmp(currBlock, "PHI")) ) {
				/* AR coefficients */
				if ( arma->p == 0 )
					if ( !(arma->AR = (sparseMatrix **) arr_alloc(maxP, sparseMatrix *)) )
						errMsg("allocation", "readARMA", "AR coefficients", 0);

				if (newBlk) {
					j = arma->p++;
					if ( !(arma->AR[j] = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
						errMsg("allocation", "readARMA", "AR[n]", 0);
					if ( !(arma->AR[j]->row = (intvec) arr_alloc(numRV*maxEps, int)) )
						errMsg("allocation", "readARMA" ,"AR[n] rows", 0);
					if ( !(arma->AR[j]->col = (intvec) arr_alloc(numRV*maxEps, int)) )
						errMsg("allocation", "readARMA" ,"AR[n] columns", 0);
					if ( !(arma->AR[j]->val = (vectorC) arr_alloc(numRV*maxEps, double)) )
						errMsg("allocation", "readARMA" ,"AR[n] coefficients", 0);
					arma->AR[j]->cnt = 0;
					newBlk = CFALSE;
				}
				arma->AR[j]->col[arma->AR[j]->cnt] = str2int(fields[0]);
				arma->AR[j]->row[arma->AR[j]->cnt] = str2int(fields[1]);
				arma->AR[j]->val[arma->AR[j]->cnt] = str2float(fields[2]);
				arma->AR[j]->cnt++;
			}
			else if ( !(strcmp(currBlock, "THETA")) ) {
				/* MA coefficients */
				if ( arma->q == 0 )
					if ( !(arma->MA = (sparseMatrix **) arr_alloc(maxQ, sparseMatrix *)) )
						errMsg("allocation", "readARMA", "MA coefficients", 0);

				if (newBlk) {
					j = arma->q++;
					if ( !(arma->MA[j] = (sparseMatrix *) mem_malloc(sizeof(sparseMatrix))) )
						errMsg("allocation", "readARMA", "MA[n]", 0);
					if ( !(arma->MA[j]->row = (intvec) arr_alloc(numRV*maxEps, int)) )
						errMsg("allocation", "readARMA" ,"MA[n] rows", 0);
					if ( !(arma->MA[j]->col = (intvec) arr_alloc(numRV*maxEps, int)) )
						errMsg("allocation", "readARMA" ,"MA[n] columns", 0);
					if ( !(arma->MA[j]->val = (vectorC) arr_alloc(numRV*maxEps, double)) )
						errMsg("allocation", "readARMA" ,"MA[n] coefficients", 0);
					arma->MA[j]->cnt = 0;
					newBlk = CFALSE;
				}
				arma->MA[j]->col[arma->MA[j]->cnt] = str2int(fields[0]);
				arma->MA[j]->row[arma->MA[j]->cnt] = str2int(fields[1]);
				arma->MA[j]->val[arma->MA[j]->cnt] = str2float(fields[2]);
				arma->MA[j]->cnt++;
			}
			else if ( !(strcmp(currBlock, "EPS")) ) {
				/* Noise parameters */
				if ( newBlk ) {
					if ( !(arma->meanEps = (vectorC) arr_alloc(maxEps, double)) )
						errMsg("allocation", "readARMA", "mean of noise process", 0);
					if ( !(arma->varEps = (vectorC) arr_alloc(maxEps, double)) )
						errMsg("allocation", "readARMA", "variance of noise process", 0);
					stoc->numGroups = 2;
					stoc->groupBeg[0] = j = 0;
					newBlk = CFALSE;
				}
				arma->meanEps[j] = str2float(fields[1]);
				arma->varEps[j] = str2float(fields[2]);
				j++; stoc->numPerGroup[1]++;
			}
			else if ( !(strcmp(currBlock, "OBS")) ) {
				/* historical observations */
				if ( newBlk ) {
					if ( !(arma->obs = (vectorC *) arr_alloc(maxT, vectorC)) )
						errMsg("allocation", "readARMA", "past observations", 0);
					arma->N = numFields;
					newBlk = CFALSE;
				}
				if ( !(arma->obs[arma->T] = (vectorC) arr_alloc(numFields, double)) )
					errMsg("allocation", "readARMA", "obs[j]", 0);
				for ( n = 0; n < numFields; n++ )
					arma->obs[arma->T][n] = str2float(fields[n]);
				arma->T++;
			}
			else if ( !(strcmp(currBlock, "NOISE")) ) {
				/* historical noise observations */
				if ( newBlk ) {
					if ( !(arma->eps = (vectorC *) arr_alloc(maxT, vectorC)) )
						errMsg("allocation", "readARMA", "past noise observations", 0);
					armaR = 0; newBlk = CFALSE;
				}
				if ( !(arma->eps[armaR] = (vectorC) arr_alloc(numFields, double)) )
					errMsg("allocation", "readARMA", "eps[j]", 0);
				for ( n = 0; n < numFields; n++ )
					arma->eps[armaR][n] = str2float(fields[n]);
				armaR++;
			}
			else if ( !(strcmp(currBlock, "TREND")) ) {
				/* time series trend */
				if ( newBlk ) {
					if ( !(arma->eta = (vectorC *) arr_alloc(maxT, vectorC)) )
						errMsg("allocation", "readARMA", "time series trend", 0);
					armaE = 0; newBlk = CFALSE;
				}
				if ( !(arma->eta[armaE] = (vectorC) arr_alloc(numFields, double)) )
					errMsg("allocation", "readARMA", "eta[j]", 0);
				for ( n = 0; n < numFields; n++ )
					arma->eta[armaE][n] = str2float(fields[n]);
				armaE++;
			}
			else if ( !(strcmp(currBlock, "SEASONALITY")) ) {
				/* time series seasonality */
				if ( newBlk ) {
					if ( !(arma->sigma = (vectorC *) arr_alloc(maxT, vectorC)) )
						errMsg("allocation", "readARMA", "time series seasonality", 0);
					armaS = 0; newBlk = CFALSE;
				}
				if ( !(arma->sigma[armaS] = (vectorC) arr_alloc(numFields, double)) )
					errMsg("allocation", "readARMA", "sigma[j]", 0);
				for ( n = 0; n < numFields; n++ )
					arma->sigma[armaS][n] = str2float(fields[n]);
				armaS++;
			}
			else {
				errMsg("read", "readARMA", "unknown field encountered", 0);
				return 1;
			}
		}
	}

	/* check the length of residual, trend and seasonality time series */
	if ( arma->T != armaR )
		printf("Warning: length of residual time series does not match with observations (%d).\n", arma->T - armaR);
	if (arma->T != armaE)
		printf("Warning: length of trend time series does not match with observations (%d).\n", arma->T - armaE);
	if ( arma->T != armaS )
		printf("Warning: length of seasonality time series does not match with observations (%d).\n", arma->T - armaS);

	stoc->arma = arma;
	/* setup mean to based on historical observation corresponding to first time period */
	n = -1;
	for ( j = 0; j < stoc->numOmega; j++ ) {
		if ( j % stoc->arma->N == 0 )
			n++;
		stoc->mean[j] = stoc->arma->obs[n][j] + stoc->arma->eta[n][j];
		if ( stoc->arma->sigma != NULL )
			stoc->mean[j] += stoc->arma->sigma[n][j];
	}

	return 0;
}//END readARMA()

int readScenarios(FILE *fptr, stringC *fields, oneProblem *orig, timeType *tim, int maxOmegas, int maxVals, stocType *stoc) {
	stringC	*rvRows, *rvCols, *scenName;
	char 	strType;
	int  	n, r, c, numFields, maxScenarios = 100, numScen = 0, parentIdx;

	/* allocate memory to hold the names of random variable */
	if ( !(rvRows = (stringC *) arr_alloc(maxOmegas, stringC)) )
		errMsg("allocation", "readScenarios", "rvNames", 0);
	if ( !(rvCols = (stringC *) arr_alloc(maxOmegas, stringC)) )
		errMsg("allocation", "readScenarios", "rvNames", 0);
	if ( !(scenName = (stringC *) arr_alloc(maxScenarios, stringC)) )
		errMsg("allocation", "readScenarios", "scenNames", 0);
	if ( !(stoc->probs[0] = (vectorC) arr_alloc(maxScenarios, double)) )
		errMsg("allocation", "readScenarios", "scenario probability", 0);

	if ( !(strcmp(fields[1], "DISCRETE")) ) {
		/* store the type of stochastic process encountered */
		sprintf(stoc->type, "SCENARIOS_DISCRETE");
		while (CTRUE) {
			NEXT_LINE: getLine(&fptr, fields, &strType, &numFields);
			if (strType != 'f')
				break; 										//Encountered ENDATA
			if ( !(strcmp(fields[0], "SC")) ) {
				/* New scenario encountered */
				if (!(strcmp(fields[2], "ROOT")) ) {
					/* The current scenario is the root scenario. Need to copy the names of random variable rows and columns. */
					if ( !(scenName[numScen] = (stringC) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readScenarios", "scenNames[n]", 0);
					strcpy(scenName[numScen], fields[1]);
					parentIdx = -1;
				}
				else {
					/* a non-root scenario encountered */
					if ( !(scenName[numScen] = (stringC) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readScenarios", "scenNames[n]", 0);
					strcpy(scenName[numScen], fields[1]);
					parentIdx = 0;
					while (parentIdx < numScen) {
						if ( !(strcmp(scenName[parentIdx], fields[2])) )
							break;
						parentIdx++;
					}
				}
				if ( parentIdx >= 0 )
					for ( n = 0; n < stoc->numOmega; n++ )
						stoc->vals[n][numScen] = stoc->vals[n][parentIdx];
				stoc->probs[0][numScen++] = str2float(fields[3]);
			}
			else {
				n = 0;
				while (n < stoc->numOmega ) {
					if ( !(strcmp(fields[0], rvCols[n])) && !(strcmp(fields[1], rvRows[n])) )
						break;
					n++;
				}
				if ( n == stoc->numOmega ) {
					/* New omega encountered. Note the row and column information. */
					/* First, begin by identifying row and column coordinates in the problem */
					if ( !(strcmp(fields[0], "RHS")) )
						c = -1;
					else {
						c = 0;
						while ( c < orig->mac ){
							if ( !(strcmp(fields[0], orig->cname[c])) )
								break;
							c++;
						}
					}
					if ( c == orig->mac ) {
						errMsg("read", "readIndep", "unknown column name in the stoch file", 0);
						return 1;
					}
					else if ( c < tim->col[1] && c != -1)
						/* variable belongs to the first stage, hence ignore */
						goto NEXT_LINE;

					if ( !(strcmp(fields[1], orig->objname)) )
						r = -1;
					else {
						r = 0;
						while (r < orig->mar ) {
							if ( !(strcmp(fields[1], orig->rname[r])) )
								break;
							r++;
						}
					}
					if ( r == orig->mar ) {
						errMsg("read", "readIndep", "unknown row name in the stoch file", 0);
						return 1;
					}
					else if ( r < tim->row[1] && r != -1)
						/* variable belongs to the first stage, hence ignore */
						goto NEXT_LINE;

					if ( !(rvRows[stoc->numOmega] = (stringC) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readIndep", "rvNames[n]", 0);
					if ( !(rvCols[stoc->numOmega] = (stringC) arr_alloc(NAMESIZE, char)) )
						errMsg("allocation", "readIndep", "rvNames[n]", 0);
					if ( !(stoc->vals[stoc->numOmega] = (vectorC) arr_alloc(maxVals, double)) )
						errMsg("allocation", "readIndep","omega.vals[n]", 0);

					strcpy(rvCols[stoc->numOmega], fields[0]);
					strcpy(rvRows[stoc->numOmega], fields[1]);
					stoc->numVals[stoc->numOmega] = 0;

					stoc->col[stoc->numOmega] = c;
					stoc->row[stoc->numOmega] = r;
					stoc->numOmega++;
				}
				stoc->vals[n][numScen-1] = str2float(fields[2]);
				stoc->numVals[n]++;
			}
		}
	}
	else {
		errMsg("read", "readScenarios", "unknown keyword in stoch file", 0);
		return 1;
	}

	/* free up memory from temporary variables */
	for ( n = 0; n < stoc->numOmega; n++) {
		if ( rvRows[n] ) mem_free(rvRows[n]);
		if ( rvCols[n] ) mem_free(rvCols[n]);
	}
	mem_free(rvRows); mem_free(rvCols);
	for ( n = 0; n < numScen; n++)
		if (scenName[n]) mem_free(scenName[n]);
	mem_free(scenName);

	/* compute the mean value for each random variable */
	for (n = 0; n < stoc->numOmega; n++ ) {
		for ( c = 0; c < numScen; c++ )
			stoc->mean[n] += stoc->probs[0][c]*stoc->vals[n][c];
	}

	/* reallocate memory to elements of stocType based on the exact sizes */
	stoc->col 		= (intvec) mem_realloc(stoc->col, stoc->numOmega*sizeof(int));
	stoc->row 		= (intvec) mem_realloc(stoc->row, stoc->numOmega*sizeof(int));
	stoc->mean      = (vectorC) mem_realloc(stoc->mean, stoc->numOmega*sizeof(double));
	stoc->numVals	= (intvec) mem_realloc(stoc->numVals, stoc->numOmega*sizeof(int));
	for ( n = 0; n < stoc->numOmega; n++ )
		stoc->vals[n] = (vectorC) mem_realloc(stoc->vals[n], numScen*sizeof(double));
	stoc->vals 		= (vectorC *) mem_realloc(stoc->vals, stoc->numOmega*sizeof(vectorC));
	stoc->probs[0] 	= (vectorC) mem_realloc(stoc->probs[0], numScen*sizeof(double));
	stoc->probs 	= (vectorC *) mem_realloc(stoc->probs, 1*sizeof(vectorC));
	mem_free(stoc->groupBeg); stoc->groupBeg = NULL;
	mem_free(stoc->numPerGroup); stoc->numPerGroup = NULL;

	return 0;
}//END readScenarios()

void freeOneProblem(oneProblem *p) {

	if(p){
		if ( p->lp ) freeProblem(p->lp);
		if(p->name) mem_free(p->name);
		if(p->objx) mem_free(p->objx);
		if(p->rhsx) mem_free(p->rhsx);
		if(p->senx) mem_free(p->senx);
		if(p->bdl) mem_free(p->bdl);
		if(p->bdu) mem_free(p->bdu);
		if(p->ctype) mem_free(p->ctype);
		if(p->matbeg) mem_free(p->matbeg);
		if(p->matval) mem_free(p->matval);
		if(p->matind) mem_free(p->matind);
		if(p->matcnt) mem_free(p->matcnt);
		if(p->objname) mem_free(p->objname);
		if(p->cname) mem_free(p->cname);
		if(p->rname) mem_free(p->rname);
		if(p->cstore) mem_free(p->cstore);
		if(p->rstore) mem_free(p->rstore);
		mem_free(p);
	}

}//END freeOneProblem()

void freeTimeType(timeType *tim) {
	int n;

	if(tim){
		if (tim->colStg) mem_free(tim->colStg);
		if (tim->rowStg) mem_free(tim->rowStg);
		if (tim->col) mem_free(tim->col);
		if (tim->row) mem_free(tim->row);
		if (tim->stgNames) {
			for (n = 0; n < tim->numStages; n++ )
				if (tim->stgNames[n]) mem_free(tim->stgNames[n]);
			mem_free(tim->stgNames);
		}
		mem_free(tim);
	}
}//END freeTime()

void freeStocType(stocType *stoc) {
	int n;

	if ( stoc ) {
		if ( stoc->col ) mem_free(stoc->col);
		if ( stoc->row ) mem_free(stoc->row);
		if ( stoc->mean ) mem_free(stoc->mean);
		if ( stoc->numVals ) mem_free(stoc->numVals);
		if ( stoc->vals) {
			if ( !(strcmp(stoc->type, "INDEP_NORMAL")) )
				mem_free(stoc->vals[0]);
			else {
				for ( n = 0; n < stoc->numOmega; n++ )
					if ( stoc->vals[n] ) mem_free(stoc->vals[n]);
			}
			mem_free(stoc->vals);
		}
		if ( stoc->probs ) {
			if ( !(strcmp(stoc->type, "SCENARIOS_DISCRETE")) ) {
				if (stoc->probs[0]) mem_free(stoc->probs[0]);
			}
			else {
				for ( n = 0; n < max(stoc->numOmega, stoc->numGroups); n++ )
					if ( stoc->probs[n] ) mem_free(stoc->probs[n]);
			}
			mem_free(stoc->probs);
		}
		if ( stoc->groupBeg) mem_free(stoc->groupBeg);
		if ( stoc->numPerGroup) mem_free(stoc->numPerGroup);
		if ( stoc->type) mem_free(stoc->type);
		if ( stoc->arma) freeARMAtype(stoc->arma);
		mem_free(stoc);
	}

}//END freeStocType()

void freeARMAtype(armaType *arma) {
	int n;

	if ( arma->AR ) {
		for ( n = 0; n < arma->p; n++ )
			if (arma->AR[n]) freeSparseMatrix(arma->AR[n]);
		mem_free(arma->AR);
	}
	if ( arma->MA ) {
		for ( n = 0; n < arma->q; n++ )
			if (arma->MA[n]) freeSparseMatrix(arma->MA[n]);
		mem_free(arma->MA);
	}
	if ( arma->obs ) {
		for ( n = 0; n < arma->T; n++ )
			if ( arma->obs[n]) mem_free(arma->obs[n]);
		mem_free(arma->obs);
	}
	if ( arma->eps ) {
		for ( n = 0; n < arma->T; n++ )
			if ( arma->eps[n]) mem_free(arma->eps[n]);
		mem_free(arma->eps);
	}
	if ( arma->eta ) {
		for ( n = 0; n < arma->T; n++ )
			if ( arma->eta[n]) mem_free(arma->eta[n]);
		mem_free(arma->eta);
	}
	if ( arma->sigma ) {
		for ( n = 0; n < arma->T; n++ )
			if ( arma->sigma[n]) mem_free(arma->sigma[n]);
		mem_free(arma->sigma);
	}
	if (arma->meanEps) mem_free(arma->meanEps);
	if (arma->varEps) mem_free(arma->varEps);

	mem_free(arma);

}//END freeARMAtype()

void writeStocType (stocType *stoc) {

	printf("-----------------------------------------------------------------------------------------------------\n");
	printf("Type of stoc file          = %s\n", stoc->type);
	printf("Number of random variables = %d\n", stoc->numOmega);
	printf("Number of groups           = %d \t with number of elements per group = ", stoc->numGroups);
	printIntvec(stoc->numPerGroup-1, stoc->numGroups, NULL);
	printf("Number of ciphers used     = %d\n", stoc->numCipher);

}//END writeStocType()
