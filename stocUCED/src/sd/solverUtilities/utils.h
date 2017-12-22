/*
 * utils.h
 *
 *  Created on: Sep 28, 2015
 *      Author: Harsha Gangammanavar
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <sys/stat.h>

#define 	NAMESIZE			32
#define		BLOCKSIZE			256
#define		MAXBITS				sizeof(int) * 8
#define		GE					'G'
#define		LE					'L'
#define		EQ					'E'

typedef 	enum {FALSE, TRUE} 	BOOL;

#define 	mem_malloc(n) 		log_alloc("malloc : " #n,malloc((n)), (n))
#define 	mem_calloc(n,size) 	log_alloc("calloc : " #n " : " #size, calloc((n),(size)), ((n) * size))
#define		arr_alloc(n,type)	(type *) mem_calloc((n),sizeof(type))
#define 	mem_realloc(ptr, n) log_realloc("realloc : " #ptr " : " #n,(ptr), realloc((ptr),(n)), (n))
#define 	mem_free(ptr) 		free(ptr)
#define 	min(X, Y) 			((X) <= (Y) ? (X) : (Y))
#define 	max(X, Y) 			((X) >= (Y) ? (X) : (Y))
#define		DBL_ABS(x)			((x) > 0.0 ? (x) : -(x))

typedef 	char				*stringC;
typedef		double				*vectorC;
typedef		int					*intvec;

typedef struct sparseVector_{
	int		cnt;
	intvec	col;
	vectorC	val;
}sparseVector;

typedef struct sparseMatrix_{
	int		cnt;
	intvec	col;
	intvec	row;
	vectorC	val;
}sparseMatrix;

/* following sub-routines can be found in utility.c */
FILE *openFile(stringC outputDir, stringC fname, char *mode);
void createOutputDir(stringC outputDir, stringC algoName, stringC Name);
void errMsg(stringC type, stringC place, stringC item, int quit);

int getLine(FILE **input, stringC *fields, char *type, int *numFields);
void copyFields(stringC *fields, int numFields, vectorC vec);
int removeSpaces (char *field);

void *log_alloc(char *string, void *return_ptr, int size);
void *log_realloc(char *string, void *free_ptr, void *alloc_ptr, int size);

double str2float(char *stringC);
int str2int(char *stringC);
int getNumBits(int num);
double oneNorm(vectorC a, int len);
double twoNorm(vectorC a, vectorC b, int len);
double calcVariance(vectorC x, int lenX);
double vXv(vectorC a, vectorC b, intvec idxCol, int len);
double vXvSparse(vectorC v, sparseVector *vSparse);
vectorC MSparsexvAdd(sparseMatrix *M, vectorC v, vectorC ans);
vectorC MSparsexvSub(sparseMatrix *M, vectorC v, vectorC ans);
vectorC vxMSparse(vectorC v, sparseMatrix *M, int len);
void vPlusv(vectorC a, vectorC b, double mult, int len);
double smooth(double newVal, double oldVal, double factor);

vectorC reduceVector(double *f_vect, int *row, int num_elem);
vectorC expandVector(vectorC red, intvec col, int redElems, int expElems);
BOOL equalVector(vectorC a, vectorC b, int len, double tolerance);
BOOL equalIntvec(intvec a, intvec b, int len);
BOOL isZeroVector(vectorC a, int len, double tolerance);
BOOL isInteger(vectorC x, int length, int startIdx, int endIdx, double tolerance);
vectorC duplicVector(double *a, int len);
intvec duplicIntvec(intvec a, int len);
void copyVector(vectorC a, vectorC b, int len, BOOL isOneNorm);
void copyIntvec (intvec a, intvec b, int len);
void addVectors(vectorC a, vectorC b, intvec indices, int len);

void trPrint(stringC routine, int type);
void printVector(vectorC vec, int len, FILE *fptr);
void printIntvec(intvec vec, int len, FILE *fptr);
void printSparseVector(vectorC vec, intvec indices, int len);
void printSparseMatrix(sparseMatrix *V, char *stringC);
void printLine();

intvec findElems(intvec allElem, int totalElem, int *numUniq);

void freeSparseMatrix(sparseMatrix *M);
void freeSparseVector(sparseVector *v);

#endif /* UTILS_H_ */
