# set the path
# set the # of scenarios
numScenarios = 2

for (rep in 0:2) {
  # set the seed
  set.seed(rep)
  
  # run the simulation
  source("runScript.R")
  
  # print scenarios
  for (i in 1:numScenarios) {
    write.csv(scenarios[,,i], sprintf("sim_rep%d_scen%d.csv", rep+1, i))
  }
}

