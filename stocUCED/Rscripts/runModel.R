# run.R
# Fits the VAR-models

# initializations
# rm(list=ls(all=TRUE))
set.seed(1)
RScriptsPath = "/Users/semihatakan/Documents/Coding Projects/Power Systems/tamingDuck/tamingDuck/stocUCED/Rscripts"
dataSetPath = "/Users/semihatakan/Documents/Coding Projects/Power Systems/tamingDuck/tamingDuck/stocUCED/datasets/3d-nrel118/stocProcesses"
simLength = 24 #days
# numScenarios = 10
lookahead = 3 #hours
  
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

### Read Data ###
dataset = readNRELData(inputDir = sprintf("%s/",dataSetPath))

### Wind Model ###
wModel = windVARModel(dataSet = dataset, decomposeBy = "monthly", identifier = "January", forecastType = "DA", lag.max = 5, infocrit = "SC")

for (name in names(wModel$model$est$varresult)) {
  if (anyNA(wModel$model$est$varresult[[name]]$coefficients)) {
    print("Error! The VAR model cannot estimate regression coefficients!")
  }
}

### Solar Model ###
sModel = solarVARModel(dataSet = dataset, decomposeBy = "monthly", identifier = "January", forecastType = "DA", lag.max = 5, infocrit = "SC", clearSky = TRUE)

for (name in names(sModel$model$est$varresult)) {
  if (anyNA(sModel$model$est$varresult[[name]]$coefficients)) {
    print("Error! The VAR model cannot estimate regression coefficients!")
  }
}

print("Modeling is completed")