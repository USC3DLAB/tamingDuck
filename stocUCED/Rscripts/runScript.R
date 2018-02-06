# Remove all except functions
# rm(list = setdiff(ls(), lsf.str()))
# cat("\014");
print("Initiating Forecasting Module");

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
simLength <- 24;
simFrequency = 4;
numScenarios <- 2;
day = 1;

set.seed(11)

# Fit model.
if (fitModel) {
	print("Fitting Forecasting Model");
	varModel <- varModelFit( paste(dataFolder, dataType, sep = ""), fileName, dataFrequency, header = TRUE, rowNames = 1, infoCrit);
}

print("Simulating Scenarios");
#Simulate scenarios
scenarios <- varSimulate(varModel, simLength, numScenarios, simFrequency, dataType)
#scenarios = scenarios[,,1]
#scenarios <- as.numeric(scenarios)
print("Forecasting is completed")
