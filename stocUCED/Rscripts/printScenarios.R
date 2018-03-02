# initializations
rm(list=ls(all=TRUE))
set.seed(0)
numScenarios = 1000
nRep = 7

# fit the model
source("runModel.R")

# print the mean scenario

# compute the mean scenario (inc. lookahead)
meanDailyWTrend <- wModel$model$dailyTrend
meanDailyWTrend <- abind::abind(meanDailyWTrend, wModel$model$dailyTrend[1:lookahead,], along = 1)
meanDailyWTrend <- array(meanDailyWTrend, c(nrow(meanDailyWTrend), ncol(meanDailyWTrend), 1)) # make it 3D array

meanDailySTrend <- sModel$model$avgSky
meanDailySTrend <- abind::abind(meanDailySTrend, sModel$model$avgSky[1:lookahead,], along = 1)
meanDailySTrend <- array(meanDailySTrend, c(nrow(meanDailySTrend), ncol(meanDailySTrend), 1)) # make it 3D array

# interpolate
meanDailyWTrend <- interpolate(simLength = simLength, lookahead = lookahead, wModel, meanDailyWTrend, 1)
meanDailySTrend <- interpolate(simLength = simLength, lookahead = lookahead, sModel, meanDailySTrend, 1)

# clean up
sunrise = (sModel$model$dayTime[1]-1-1)*4+2
sunset = sModel$model$dayTime[ length(sModel$model$dayTime) ]*4
meanDailySTrend[-(sunrise:sunset),,] = 0

meanDailyWTrend <- meanDailyWTrend[1:96,,1]
meanDailySTrend <- meanDailySTrend[1:96,,1]

#replicate
nrep = dim(wModel$ts)[3]

meanDailyWTrend = matrix(rep(t(meanDailyWTrend),nrep), ncol=ncol(meanDailyWTrend), byrow=TRUE)
meanDailySTrend = matrix(rep(t(meanDailySTrend),nrep), ncol=ncol(meanDailySTrend), byrow=TRUE)

# print
write.csv(meanDailyWTrend, "DA_WIND.csv")
write.csv(meanDailySTrend, "DA_SOLAR.csv")

# simulate
for (rep in 1:nRep) {
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

