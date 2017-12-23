/*
 * readFIles.c
 *
 *  Created on: Oct 17, 2016
 *      Author: Harsha Gangammanavar
 *              Southern Methodist University.
 *  	Report: harsha@smu.edu
 */

#include "utils.h"
#include "solver.h"
#include "smps.h"
#include "prob.h"

// global variables
long int MEM_USED;
string 	outputDir;

// subroutines used
void parseCmdLine(int argc, string *argv, string algoName, string probName, string inputDir);
void setupDir(string algoName, string probName);

int main(int argc, char *argv[]) {
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;
	probType **prob = NULL;
	vector	meanSol, lb;
	char probName[NAMESIZE], algoName[NAMESIZE], inputDir[2*BLOCKSIZE];

	/* parse command line to obtain input from user regarding problem name and the algorithm for which the problem is being read for. */
	parseCmdLine(argc, argv, algoName, probName, inputDir);

	/* set up input and output directories for algorithm and problem being solved */
	sprintf(inputDir, "%s/", inputDir);
	outputDir = (string) mem_malloc(BLOCKSIZE*sizeof(char));
	setupDir(algoName, probName);

	/* open solver environment */
	openSolver();

	/* read problem information from the SMPS files */
	if ( readFiles(inputDir, probName, &orig, &tim, &stoc) ) {
		errMsg("read", "main", "failed to read problem files", 0);
		goto TERMINATE;
	}
	printf("Successfully read '%s' SMPS files for '%s' algorithm.\n", probName, algoName);

#ifdef INPUT_CHECK
	/* print the summary of problems that have been read */
	printReadSummary(orig, tim, stoc);
#endif

	/* setup mean value problem which will act as reference for all future computations */
	meanSol = meanProblem(orig, stoc);
	if ( meanSol == NULL ) {
		errMsg("setup", "setupAlgo", "failed to setup and solve mean value problem", 0);
		return 1;
	}

	/* calculate lower bounds for each stage */
	lb = calcLowerBound(orig, tim, stoc);
	if ( lb == NULL ) {
		errMsg("setup", "setupAlgo", "failed to compute lower bounds on stage problem", 0);
		return 1;
	}

	/* decompose the problem into master and subproblem */
	prob = newProb(orig, stoc, tim, lb, 0.0001);
	if ( prob == NULL ) {
		errMsg("setup", "setupAlgo", "failed to update probType with elements specific to algorithm", 0);
		goto TERMINATE;
	}
	printf("Successfully decomposed the problem '%s' for '%s' algorithm.\n", probName, algoName);

#ifdef DECOMPOSE_CHECK
	printDecomposeSummary(tim, prob);
#endif

	TERMINATE:
	/* release problem structures */
	freeProbType(prob, tim->numStages);
	freeOneProblem(orig);
	freeStocType(stoc);
	freeTimeType(tim);
	mem_free(outputDir); mem_free(meanSol); mem_free(lb);
	/* close solver environment and release all structures */
	closeSolver();

	return 0;
}//END

/* Parse the command line to obtain the names of the algorithm and the problem */
void parseCmdLine(int argc, string *argv, string algoName, string probName, string inputDir) {

	switch (argc) {
	case 4:
		strcpy(algoName, argv[1]);
		strcpy(probName, argv[2]);
		strcpy(inputDir, argv[3]);
		break;
	case 3:
		strcpy(algoName, argv[1]);
		strcpy(probName, argv[2]);
		printf("Please enter the input directory: ");
		scanf("%s", inputDir);
		break;
	case 2:
		strcpy(algoName, argv[1]);
		printf("Please enter the name of the problem: ");
		scanf("%s", probName);
		printf("Please enter the input directory: ");
		scanf("%s", inputDir);
		break;
	default:
		printf("Please enter the name of the algorithm: ");
		scanf("%s", algoName);
		printf("Please enter the name of the problem: ");
		scanf("%s", probName);
		printf("Please enter the input directory: ");
		scanf("%s", inputDir);
	}

}//END parseCmdLine

/* setup an output directory for the problem in the algorithms directory */
void setupDir(string algoName, string probName) {
	char buffer[2*BLOCKSIZE];

	sprintf(outputDir, "../../spOutput/%s/", algoName);
	sprintf(buffer, "mkdir %s", outputDir);
	system(buffer);
	sprintf(outputDir, "../../spOutput/%s/%s/", algoName, probName);
	sprintf(buffer, "mkdir %s", outputDir);
	if (system(buffer)){
		sprintf(buffer, "rm -r %s*", outputDir);
		system(buffer);
	}

}//END setuDir()
