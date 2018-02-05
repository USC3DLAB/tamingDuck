# Remove all except functions
rm(list = setdiff(ls(), lsf.str()))
cat("\014");
options(warn=-1);

# Data parameters
dataFolder = "/home/gjharsha/Documents/workspace/tamingDuck/stocUCED/datasets/3d-nrel118/Solar"
fileName   = "DA_1.csv";
dataFrequency = 1;
infoCrit = "SC"
type = "Solar"

# Simulation parameters 
simLength <- 24;
simFrequency = 4;
numScenarios <- 1;
day = 1;

set.seed(11)

# Fit model.
varModel <- varModelFit(dataFolder, fileName, dataFrequency, header = TRUE, rowNames = 1, infoCrit);

#Simulate scenarios
scenarios <- varSimulate(varModel, simLength, numScenarios, simFrequency, type)
