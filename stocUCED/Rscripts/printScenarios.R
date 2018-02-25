# initializations
rm(list=ls(all=TRUE))
set.seed(0)
numScenarios = 1000
nRep = 7

# fit the model
source("runModel.R")

# simulate
for (rep in 5:nRep) {
  # set the seed
  set.seed(rep)
  
  # run the simulation
  source("runSimulate.R")

  # combine sample paths
  scenarios = abind::abind(sPaths, wPaths, along=2)
  
  # print scenarios
  for (i in 1:numScenarios) {
    write.csv(scenarios[,,i], sprintf("sim_rep%d_scen%d.csv", rep, i))
  }
  
  # rm scenarios
  rm(wPaths)
  rm(sPaths)
  rm(scenarios)
}

