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

typedef 	char				*string;
typedef		double				*vector;
typedef		int					*intvec;
typedef 	enum {FALSE, TRUE} 	BOOL;

#define 	mem_malloc(n) 		log_alloc("malloc : " #n,malloc((n)), (n))
#define 	mem_calloc(n,size) 	log_alloc("calloc : " #n " : " #size, calloc((n),(size)), ((n) * size))
#define		arr_alloc(n,type)	(type *) mem_calloc((n),sizeof(type))
#define 	mem_realloc(ptr, n) log_realloc("realloc : " #ptr " : " #n,(ptr), realloc((ptr),(n)), (n))
#define 	mem_free(ptr) 		free(ptr)
#define 	min(X, Y) 			((X) <= (Y) ? (X) : (Y))
#define 	max(X, Y) 			((X) >= (Y) ? (X) : (Y))
#define		DBL_ABS(x)			((x) > 0.0 ? (x) : -(x))

typedef struct sparseVector_{
	int		cnt;
	intvec	col;
	vector	val;
}sparseVector;

typedef struct sparseMatrix_{
	int		cnt;
	intvec	col;
	intvec	row;
	vector	val;
}sparseMatrix;

/* following sub-routines can be found in utility.c */
FILE *openFile(string outputDir, string fname, char *mode);
void createOutputDir(string outputDir, string algoName, string probName);
void errMsg(string type, string place, string item, int quit);

int getLine(FILE **input, string *fields, char *type, int *numFields);
void copyFields(string *fields, int numFields, vector vec);
int removeSpaces (char *field);

void *log_alloc(char *string, void *return_ptr, int size);
void *log_realloc(char *string, void *free_ptr, void *alloc_ptr, int size);

double str2float(char *string);
int str2int(char *string);
int getNumBits(int num);
double oneNorm(vector a, int len);
double twoNorm(vector a, vector b, int len);
double calcVariance(vector x, int lenX);
double vXv(vector a, vector b, intvec idxCol, int len);
double vXvSparse(vector v, sparseVector *vSparse);
vector MSparsexvAdd(sparseMatrix *M, vector v, vector ans);
vector MSparsexvSub(sparseMatrix *M, vector v, vector ans);
vector vxMSparse(vector v, sparseMatrix *M, int len);
void vPlusv(vector a, vector b, double mult, int len);
double smooth(double new, double old, double factor);

vector reduceVector(double *f_vect, int *row, int num_elem);
vector expandVector(vector red, intvec col, int redElems, int expElems);
BOOL equalVector(vector a, vector b, int len, double tolerance);
BOOL equalIntvec(intvec a, intvec b, int len);
BOOL isZeroVector(vector a, int len, double tolerance);
BOOL isInteger(vector x, int length, int startIdx, int endIdx, double tolerance);
vector duplicVector(double *a, int len);
intvec duplicIntvec(intvec a, int len);
void copyVector(vector a, vector b, int len, BOOL isOneNorm);
void copyIntvec (intvec a, intvec b, int len);
void addVectors(vector a, vector b, intvec indices, int len);

void trPrint(string routine, int type);
void printVector(vector vec, int len, FILE *fptr);
void printIntvec(intvec vec, int len, FILE *fptr);
void printSparseVector(vector vec, intvec indices, int len);
void printSparseMatrix(sparseMatrix *V, char *string);
void printLine();

intvec findElems(intvec allElem, int totalElem, int *numUniq);

void freeSparseMatrix(sparseMatrix *M);
void freeSparseVector(sparseVector *v);

#endif /* UTILS_H_ */
