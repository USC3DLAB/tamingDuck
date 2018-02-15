/*
 * utility.c
 *
 *  Created on: Apr 20, 2014
 *      Author: Harsha Gangammanavar
 */

#include "utils.h"

extern long MEM_USED;

FILE *openFile(stringC outputDir, stringC fname, char *mode) {
	FILE *fptr;
	char buffer[2*BLOCKSIZE];

	strcpy(buffer, outputDir);
	strcat(buffer, fname);

	fptr = fopen(buffer, mode);
	if ( fptr == NULL ) {
		fprintf(stderr, "failed to open file %s", fname);
		return NULL;
	}

	return fptr;
}//END openFile()

void createOutputDir(stringC outputDir, stringC algoName, stringC probName) {
	struct stat st;
	char buffer[2*BLOCKSIZE];

	strcat(outputDir,algoName);
	strcat(outputDir,"/");
	if ( stat(outputDir, &st) ) {
		sprintf(buffer, "mkdir %s", outputDir);
		system(buffer);
	}
	strcat(outputDir, probName);
	strcat(outputDir, "/");
	if ( stat(outputDir, &st) ) {
		sprintf(buffer, "mkdir %s", outputDir);
		system(buffer);
	}
	else {
		sprintf(buffer, "rm -r %s*", outputDir);
		system(buffer);
	}

}//END createOutputDir()

void errMsg(stringC type, stringC place, stringC item, int quit){
	fprintf(stderr, "\nError :: Type - %s;  Function - %s(); Item - %s\n", type, place, item);
	if (quit)
		exit(1);
}//err_msg()

/* The function getLine() reads an input line, checks its length, and determines the appropriate parsing function to call. In the case of a line consisting of only a single carriage return, the function
 * continues to read in lines until a stringC of non_zero length or an EOF is encountered. Upon return of fields from the parsing functions, each field will be sent to remove_spaces for the removal of all
 * blank spaces.  The getLine function then returns control to the calling function. */
int getLine(FILE **input, stringC *fields, char *type, int *numFields) {
	char 	input_str[BLOCKSIZE], *strptr, *token;
	long	len = 0, n, stat;

	input_str[0] = '*';

	/*  test for and skip over empty lines and comments  */
	while (len <= 1 || input_str[0] == '*') {
		strptr = fgets(input_str, 10000, *input);
		if (strptr == NULL)
			return 1;
		len = strlen(input_str);
	}

	/*  identify type of input stringC (title or field)  */
	if (input_str[0] >= '0' && input_str[0] <= 'Z') {
		/* input stringC is title stringC, marked by type = 't' */
		sscanf(input_str, "%s %s", fields[0], fields[1]);
		type[0] = 't';
		n = 2;
	}
	else {
		/* input stringC is field stringC, marked by type = 'f' */
		token = strtok(input_str, " \t");
		n = 0;
		while( token != NULL ) {
			strcpy(fields[n++], token);
			token = strtok(NULL, " \t");
		}
		type[0] = 'f';
	}

	(*numFields) = n;
	for ( n = 0; n < (*numFields); n++ ) {
		stat = removeSpaces(fields[n]);
		if ( stat == 0 ) {
			(*numFields)--;
		}
	}

	return 0;
}//END getLine()

/* The function removeSpace() removes additional spaces from a	stringC. This function will be used after an input stringC has been broken into its appropriate fields according to column values. */
int removeSpaces (char *field) {
	char *p, *q, len = 0;;

	p = field;
	q = p;
	while (*p != '\0' && *p != '\n') {
		if (*p >= 33) {
			*q++ = *p++;
			len++;
		}
		else
			p++;
	}
	*q = '\0';

	return len;
}//END removeSpaces

void trPrint(stringC routine, int type){
	if ( type == 1 )
		printf("Entering :: %s()\n", routine);
	else
		printf("Exiting  :: %s()\n", routine);
}//END trPrint()

void *log_alloc(char *stringC, void *return_ptr, int size) {
	MEM_USED += size;
	return return_ptr;
}//END log_alloc()


void *log_realloc(char *stringC, void *free_ptr, void *alloc_ptr, int size) {
	MEM_USED += size;
	return alloc_ptr;
}//END log_realloc()


double str2float(char *stringC){
	double val;
	sscanf(stringC, "%lf", &val);
	return val;
}//END str_to_float()

int str2int(char *stringC) {
	int val;
	sscanf(stringC, "%d", &val);
	return val;
}//END str2int()

/* This function returns the minimum number of bits needed to represent a given number. */
int getNumBits(int num) {
  int 	hi_bit = 1;
  int 	numBits;

  for (numBits = 0; hi_bit <= num; numBits++)
    hi_bit = hi_bit << 1;

  return numBits;
}//END getNumBits()


double oneNorm(vectorC a, int len) {
	int		cnt;
	double	sum;

	sum = 0.0;
	for (cnt = 0; cnt < len; cnt++)
		sum += DBL_ABS(a[cnt]);

	return sum;
}//END oneNorm()

double twoNorm(vectorC a, vectorC b, int len) {
	int 	cnt;
	double	norm = 0.0;

	if (b != NULL)
		for (cnt = 1; cnt <= len; cnt++ )
			norm += pow((a[cnt]-b[cnt]), 2);
	else

	norm = sqrt(norm);
	return norm;
}//END twoNorm()

double calcVariance(vectorC x, int lenX) {
    double 	mean, vari, temp;
    int 	cnt;

    temp = 0.0;
    vari = 0.0; mean = x[0];
    for (cnt = 1; cnt < lenX; cnt++) {
        temp = mean;
        mean = mean + (x[cnt] - mean) / (double) (cnt + 1);
        vari = (1 - 1 / (double) cnt) * vari + (cnt + 1) * (mean - temp) * (mean - temp);
    }

    return vari;
}//END calcVariance()

double vXv(vectorC a, vectorC b, intvec idxCol, int len) {
	double ans = 0.0;
	int n;

	if(idxCol == NULL)
		for ( n = 1; n <= len; n++ )
			ans += a[n]*b[n];
	else
		for (n = 1; n <= len; n++ )
			ans += a[n]*b[idxCol[n]];

	return ans;
}//END vXv()

double vXvSparse(vectorC v, sparseVector *vSparse){
	int		cnt;
	double 	ans;

	ans = 0.0;
	for (cnt = 1; cnt <= vSparse->cnt; cnt++)
		ans += vSparse->val[cnt] * v[vSparse->col[cnt]];

	return ans;
}//END vXvSparse()

vectorC MSparsexvAdd(sparseMatrix *M, vectorC v, vectorC ans){
	int	cnt;

	for (cnt = 1; cnt <= M->cnt; cnt++)
		ans[M->row[cnt]] += M->val[cnt] * v[M->col[cnt]];

	return ans;
}//END MSparsexv()

vectorC MSparsexvSub(sparseMatrix *M, vectorC v, vectorC ans){
	int	cnt;

	for (cnt = 1; cnt <= M->cnt; cnt++)
		ans[M->row[cnt]] -= M->val[cnt] * v[M->col[cnt]];

	return ans;
}//END MSparsexvSub()

vectorC vxMSparse(vectorC v, sparseMatrix *M, int len) {
	int		cnt;
	vectorC	ans;

	if(!(ans = (vectorC) arr_alloc(len+1, double)))
		errMsg("allocation", "vxMSparse", "ans", 1);

	for (cnt = 1; cnt <= M->cnt; cnt++)
		ans[M->col[cnt]] += v[M->row[cnt]] * M->val[cnt];
	ans[0] = oneNorm(ans+1, len);

	return ans;
}//END PIxT()

void vPlusv(vectorC a, vectorC b, double mult, int len){
	int 	cnt;

	for ( cnt = 1; cnt <= len; cnt++ )
		a[cnt] += mult*b[cnt];
	a[0] = oneNorm(a+1, len);

}//END vPlusv()

double smooth(double new_val, double old, double factor) {
	return factor*new_val + (1-factor)*old;
}//END smooth();

vectorC reduceVector(vectorC f_vect, int *row, int num_elem){
	int		cnt;
	double 	*s_vect;

	if(!(s_vect = (vectorC) arr_alloc(num_elem+1, double)))
		errMsg("allocation", "reduceVector", "s_vect", 1);

	for (cnt = 1; cnt <= num_elem; cnt++)
		s_vect[cnt] = f_vect[row[cnt]];
	s_vect[0] = oneNorm(s_vect+1, num_elem);

	return s_vect;
}//END reduceVector()

vectorC expandVector(vectorC red, intvec col, int redElems, int expElems){
	int 	n;
	vectorC 	exp;

	if (!(exp = (vectorC) arr_alloc(expElems+1, double)) )
		errMsg("allocation", "expandVector", "expanded vectorC", 0);

	for (n = 1; n <= redElems; n++ )
		exp[col[n]] = red[n];
	exp[0] = oneNorm(exp+1, expElems);

	return exp;
}//END expandVector

BOOL equalVector(vectorC a, vectorC b, int len, double tolerance) {
	int		cnt;

	for (cnt = 1; cnt <= len; cnt++)
		if ( DBL_ABS(a[cnt] - b[cnt]) > tolerance )
			return CFALSE;
    
	return CTRUE;
}//END equalVector()

BOOL equalIntvec(intvec a, intvec b, int len) {
	int		cnt;

	for (cnt = 1; cnt <= len; cnt++)
		if ( a[cnt] != b[cnt] )
			return CFALSE;

	return CTRUE;
}//END equalIntvec()

BOOL isZeroVector(vectorC a, int len, double tolerance) {
	int		cnt;

	for (cnt = 0; cnt < len; cnt++) {
		if ( DBL_ABS(a[cnt]) >= tolerance )
			return CFALSE;
//		else
//			a[cnt] = 0.0;
	}

	return CTRUE;
}//END equalVector()

/*This function will check if a vectorC is integer with a predefined gap */
BOOL isInteger(vectorC x, int length, int startIdx, int endIdx, double tolerance){
	int i;

	for (i = startIdx+1; i < endIdx; i++)
		if (fabs(x[i] - round(x[i])) > tolerance)
			return CFALSE;

	return CTRUE;
}//END isInteger()


vectorC duplicVector(vectorC a, int len) {
	int		i;
	vectorC	b;

	if ((b = (vectorC) arr_alloc(len+1, double))) {
		for (i = 1; i <= len; i++)
			b[i] = a[i];
		b[0] = oneNorm(b+1, len);
	}
	else
		errMsg("allocation", "duplicArray", "b", 1);

	return b;
}//END duplicArray()

intvec duplicIntvec(intvec a, int len) {
	int		i;
	intvec	b;

	if ((b = (intvec) arr_alloc(len+1, int))) {
		for (i = 1; i <= len; i++)
			b[i] = a[i];
	}
	else
		errMsg("allocation", "duplicArray", "b", 1);

	return b;
}//END duplicArray()

void copyVector(vectorC a, vectorC b, int len, BOOL isOneNorm){
	int n;

	if (isOneNorm)
		for ( n = 0; n <= len; n++ )
			b[n] = a[n];
	else {
		for ( n = 1; n <= len; n++ )
			b[n] = a[n-1];
		b[0] = oneNorm(b+1, len);
	}

}//END copyVector()

void copyIntvec (intvec a, intvec b, int len) {
	int n;

	for ( n = 0; n < len; n++ )
		b[n] = a[n];

}//END copyVector()

void addVectors(vectorC a, vectorC b, intvec indices, int len){
	int n;

	if ( indices == NULL ) {
		for ( n = 1; n <= len; n++ )
			a[n] += b[n];
	}
	else {
		for ( n = 1; n <= len; n++ )
			a[indices[n]] += b[n];
	}
	a[0] = oneNorm(a+1, len);

}//END copy_arr()

void printVector(vectorC vec, int len, FILE *fptr){
	int n;

	if ( fptr == NULL ) {
		for ( n = 1; n <= len; n++ )
			printf("%4.3lf ", vec[n]);
		printf("\n");
	}
	else {
		for ( n = 1; n <= len; n++ )
			fprintf(fptr, "%4.3lf\t", vec[n]);
		fprintf(fptr, "\n");
	}

}//END printVector()

void printVectorWName(vectorC vec, stringC *vecName, int len, FILE *fptr){
	int n;

	for ( n = 1; n <= len; n++ ) {
		fprintf(fptr, "%s\t\t%4.3lf\n ", vecName[n-1],vec[n]);
		fprintf(fptr, "\n");
	}

}//END printVectorWName()

void printIntvec(intvec vec, int len, FILE *fptr){
	int n;

	if (fptr == NULL) {
		for ( n = 1; n <= len; n++ )
			printf("%d ", vec[n]);
		printf("\n");
	}
	else {
		for ( n = 1; n <= len; n++ )
			fprintf(fptr, "%d ", vec[n]);
		fprintf(fptr, "\n");
	}

}//END printIntvec()

void printSparseVector(vectorC vec, intvec indices, int len) {
	int n;

	for ( n = 1; n <= len; n++ )
		printf("%4.3lf", vec[indices[n]]);
	printf("\n");

}//END printSparseVector()

void printSparseMatrix(sparseMatrix *V, char *stringC) {
	int 	cnt;
	printf("%s (%d) ::\n\t\n", stringC, V->cnt);
	for (cnt = 1; cnt <= V->cnt; cnt++){
		printf("(%d, %d, %.2f)\t\n", V->row[cnt], V->col[cnt], V->val[cnt]);
		if ( cnt % 5 == 0 )
			printf("\n");
	}
	printf("\n");
}// END print_sparseMatrix()

void printLine() {

    printf("-------------------------------------------------------------------------- \n");

}//END printLine

intvec findElems(intvec allElem, int totalElem, int *numUniq){
	intvec	elemUniq;
	int		n, m, len;

	if(!(elemUniq = arr_alloc(totalElem+1, int)))
		errMsg("allocation", "find_cols", "colUniq", 1);

	len = 0;
	/* Copy over all the distinct non-zero elements of allElem */
	for ( n = 1; n <= totalElem; n++ ) {
		if ( allElem[n] > 0 ) {
			m = 1;
			while ( m <= len ) {
				if ( allElem[n] == elemUniq[m] )
					break;
				m++;
			}
			if ( m == len+1 )
				elemUniq[++len] = allElem[n];
		}
	}

	/* Shrink the array down to the number of distinct elements found */
	elemUniq = (intvec) mem_realloc(elemUniq, (len+1)*sizeof(int));

	elemUniq[0] = 0;
	*numUniq = len;

	return elemUniq;
}//END fin_cols

void freeSparseMatrix(sparseMatrix *M) {

	if (M->col) mem_free(M->col);
	if (M->row) mem_free(M->row);
	if (M->val) mem_free(M->val);
	mem_free(M);

}//END freeSparseMatrix()

void freeSparseVector(sparseVector *v) {

	if (v->col) mem_free(v->col);
	if (v->val) mem_free(v->val);
	mem_free(v);

}//END freeSparseMatrix()
