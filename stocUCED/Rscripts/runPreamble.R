# run.R
# Fits the VAR-models

# initializations
#rm(list=ls(all=TRUE))
#set.seed(0)
RScriptsPath = "/Users/semihatakan/Documents/Coding Projects/Power Systems/tamingDuck/tamingDuck/stocUCED/Rscripts"
dataSetPath = "/Users/semihatakan/Documents/Coding Projects/Power Systems/tamingDuck/tamingDuck/stocUCED/datasets/3d-nrel118/stocProcesses"
simLength = 24 #days
lookahead = 3 #hours

# R output re-direction
sink("R.log", append = FALSE, split = TRUE) # removes the previous log
#sink("R.log", append = TRUE, split = TRUE)

# import libraries
library("abind");

# source functions
source(sprintf("%s/readers/readNRELData.R", RScriptsPath))
source(sprintf("%s/solarModels/solarVarModel.R", RScriptsPath))
source(sprintf("%s/solarModels/solarSimulate.R", RScriptsPath))
source(sprintf("%s/windModels/windVarModel.R", RScriptsPath))
source(sprintf("%s/windModels/windSimulate.R", RScriptsPath))
source(sprintf("%s/tsUtilities/decomposeTSData.R", RScriptsPath))
source(sprintf("%s/tsUtilities/varSimulate.R", RScriptsPath))
source(sprintf("%s/tsUtilities/interpolate.R", RScriptsPath))
source(sprintf("%s/tsUtilities/prepareTSData.R", RScriptsPath))
source(sprintf("%s/tsUtilities/extractTSData.R", RScriptsPath))

### Read Data ###
dataset = readNRELData(inputDir = sprintf("%s/",dataSetPath))

print("Data is read")