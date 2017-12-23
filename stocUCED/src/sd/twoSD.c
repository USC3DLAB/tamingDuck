/*
 * twoSD.c
 *
 *  Created on: Jul 6, 2017
 *      Author: Harsha Gangammanavar
 * Institution: Southern Methodist University
 *  
 * Please send you comments or bug report to harsha (at) smu (dot) edu
 *
 */


#include <twoSD.h>

long long	MEM_USED = 0;	/* Amount of memory allocated each iteration */
stringC   	outputDir;		/* output directory */
configType	config;			/* algorithm tuning parameters */

int main (int argc, char *argv[]) {
	int 	status;
	char 	inputDir[2*BLOCKSIZE], probName[NAMESIZE];
	oneProblem *orig = NULL;
	timeType *tim = NULL;
	stocType *stoc = NULL;

	/* open solver environment */
	openSolver();

	/* read algorithm configuration files */
	status = readConfig(inputDir);
	if ( status ) {
		errMsg("read", "main", "failed to read algorithm configuration file", 0);
		goto TERMINATE;
	}

	/* read problem information */
	/* request for problem name to be solved, the path should be provided in the configuration file */
	if ( argc < 2 )
		parseCmdLine(probName);
	else
		strcpy(probName, argv[1]);

	/* read problem SMPS input files */
	status = readFiles(inputDir, probName, &orig, &tim, &stoc);
	if ( status ) {
		errMsg("read", "main", "failed to read problem files using SMPS reader", 0);
		goto TERMINATE;
	}

	/* set up output directory: using the outputDir in config file and the input problem name */
	createOutputDir(outputDir, "twoSD", probName);

	/* launch the algorithm */
	status = algo(orig, tim, stoc, inputDir, probName);
	if ( status ) {
		errMsg("allocation", "main", "failed to solve the problem using SDDP", 0);
		goto TERMINATE;
	}

	/* release structures and close solver environment */
	TERMINATE:
	freeOneProblem(orig);
	freeTimeType(tim);
	freeStocType(stoc);
	closeSolver();

	return 0;
}//END main()

void parseCmdLine(stringC probName) {

	printf("Please enter the name of the problem: ");
	scanf("%s", probName);

}//END parseCmdLine

int readConfig(stringC inputDir) {
	FILE 	*fptr;
	char	line[2*BLOCKSIZE], comment[2*BLOCKSIZE];
	int 	status;

	fptr = fopen("config.sd", "r");
	if ( fptr == NULL ) {
		errMsg("read", "readConfig", "failed to open configuration file", 0);
		return 1;
	}

	if ( !(outputDir = (stringC) mem_malloc(BLOCKSIZE*sizeof(char))) )
		errMsg("allocation", "readConfig", "outputDir", 0);

	while ((status = (fscanf(fptr, "%s", line) != EOF))) {
		if (!(strcmp(line, "INPUTDIR")))
			fscanf(fptr, "%s", inputDir);
		else if (!(strcmp(line, "OUTPUTDIR")))
			fscanf(fptr, "%s", outputDir);
		else if (!(strcmp(line, "RUN_SEED")))
			fscanf(fptr, "%lld", &config.RUN_SEED);
		else if (!(strcmp(line, "TOLERANCE")))
			fscanf(fptr, "%lf", &config.TOLERANCE);
		else if (!(strcmp(line, "MIN_ITER")))
			fscanf(fptr, "%d", &config.MIN_ITER);
		else if (!(strcmp(line, "MAX_ITER")))
			fscanf(fptr, "%d", &config.MAX_ITER);
		else if (!(strcmp(line, "MASTERTYPE")))
			fscanf(fptr, "%d", &config.MASTERTYPE);
		else if (!(strcmp(line, "CUT_MULT")))
			fscanf(fptr, "%d", &config.CUT_MULT);
		else if (!(strcmp(line, "TAU")))
			fscanf(fptr, "%d", &config.TAU);
		else if (!(strcmp(line, "MIN_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MIN_QUAD_SCALAR);
		else if (!(strcmp(line, "MAX_QUAD_SCALAR")))
			fscanf(fptr, "%lf", &config.MAX_QUAD_SCALAR);
		else if (!(strcmp(line, "R1")))
			fscanf(fptr, "%lf", &config.R1);
		else if (!(strcmp(line, "R2")))
			fscanf(fptr, "%lf", &config.R2);
		else if (!(strcmp(line, "R3")))
			fscanf(fptr, "%lf", &config.R3);
		else if (!(strcmp(line, "PI_EVAL_START")))
			fscanf(fptr, "%d", &config.PI_EVAL_START);
		else if (!(strcmp(line, "PI_CYCLE")))
			fscanf(fptr, "%d", &config.PI_CYCLE);
		else if (!(strcmp(line, "SCAN_LEN")))
			fscanf(fptr, "%d", &config.SCAN_LEN);
		else if (!(strcmp(line, "EVAL_FLAG")))
			fscanf(fptr, "%d", &config.EVAL_FLAG);
		else if (!(strcmp(line, "EVAL_SEED")))
			fscanf(fptr, "%lld", &config.EVAL_SEED);
		else if (!(strcmp(line, "EVAL_MIN_ITER")))
			fscanf(fptr, "%d", &config.EVAL_MIN_ITER);
		else if (!(strcmp(line, "EVAL_ERROR")))
			fscanf(fptr, "%lf", &config.EVAL_ERROR);
		else if (!(strcmp(line, "PRE_EPSILON")))
			fscanf(fptr, "%lf", &config.PRE_EPSILON);
		else if (!(strcmp(line, "EPSILON")))
			fscanf(fptr, "%lf", &config.EPSILON);
		else if (!(strcmp(line, "BOOTSTRAP_REP")))
			fscanf(fptr, "%d", &config.BOOTSTRAP_REP);
		else if (!strcmp(line, "//"))
			fgets(comment, 2*BLOCKSIZE, fptr);
		else {
			printf ("%s\n", line);
			errMsg("read", "readConfig", "unrecognized parameter in configuration file", 1);
		}
	}

	fclose(fptr);

	return 0;
}//END readConfig()
