# Remove all except functions
# rm(list = setdiff(ls(), lsf.str()))
# cat("\014");
print("Initiating simulation scripts");

# load packages
options(warn=-1);
suppressPackageStartupMessages(library("zoo"));		# just to shut R up
suppressPackageStartupMessages(library("sandwich"));
suppressPackageStartupMessages(library("MTS"));
suppressPackageStartupMessages(library("vars"));
library("zoo");
library("strucchange");
library("sandwich");
library("MASS");
library("urca");
library("lmtest");
library("vars");
library("MTS");
library("abind");

# Source the functions
source("preprocessNRELdata.R");
source("varModelFit.R");
source("varSimulate.R");

# Data parameters
# dataFolder = "/Users/semihatakan/Documents/Coding Projects/Power Systems/tamingDuck/tamingDuck/stocUCED/datasets/3d-nrel118/"
# dataType   = "Solar";
# fileName   = "DA.csv";
# fitModel	 = TRUE;
dataFrequency = 1;
infoCrit = "SC"

# Simulation parameters 
# simLength <- 24;
simFrequency = 4;
# numScenarios <- 2;
day = 1;

set.seed(11)

# Fit model
if (fitModel) {
	varModel = NULL;
	colNames = NULL;
	for (i in 1:length(dataTypes)) {
		print(sprintf("Fitting time-series model for %s", dataTypes[i]));
		varModel[[i]] <- varModelFit( paste(dataFolder, dataTypes[i], sep = ""), fileName, dataFrequency, header = TRUE, rowNames = 1, infoCrit);
		colNames = append(colNames, colnames(varModel[[i]]$ts$data));
	}
}

#Simulate scenarios
print("Simulating Scenarios")
numComponents = 0;
for (i in 1:length(dataTypes)) {
  numComponents = numComponents + varModel[[i]]$ts$N;
}

scenarios = array(0, dim = c(simLength*simFrequency, numComponents, numScenarios))
#scenarios = vector(mode = "double", length = numComponents*simLength*simFrequency*numScenarios)

j = 1;
for (i in 1:length(dataTypes)) {
  tmp <- varSimulate(varModel[[i]], simLength, numScenarios, simFrequency, dataTypes[i])
  # tmp[time, generator, scenario]

  scenarios[, j:(j+varModel[[i]]$ts$N-1), ] = tmp;
  j = j + varModel[[i]]$ts$N;
}
print("Simulation is completed")
