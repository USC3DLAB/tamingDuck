# initializations
rm(list=ls(all=TRUE))
set.seed(0)
numScenarios = 1000
nRep = 7
today = 45  # Feb 14, we will plan, starting with this day...  
#today = 197 # Jul 15"
forecastWindow = 14 # we will look 14 days of history for fitting the models...

source("runPreamble.R")

# print deterministic elements of the time-series
printTS <- function(srcType, fileType) {
  nReplications = nRep+1  # give it another day for look-ahead
  data = prepareTSData(dataset, srcType, "extract", fileType, today = today+nReplications, nReplications)  
  data = apply(data$ts, 2, c)
  data = interpolate(simLength = dim(data)[1], lookahead = 0, dim(data)[2], data, 1)
  if (srcType == "Solar") { # remove the night-time generations
    hourlyData = prepareTSData(dataset, srcType, "extract", fileType, today = today+nReplications, nReplications);
    hourlyData = hourlyData$ts
    for (g in 1:dim(data)[2]) {
      indices <- NULL
      for (d in 1:nReplications) {
        dayTime = which(hourlyData[,g,d] > 0)
        sunrise = (dayTime[1]-1-1)*4+2
        sunset = dayTime[ length(dayTime) ]*4
        indices = abind(indices, (sunrise:sunset) + (d-1)*96)
      }
      data[-indices,g,] = 0
    }
  }
  write.csv(data[,,1], sprintf("%s_%s.csv", fileType, srcType))
}
printTS("Load", "DA")
printTS("Load", "RT")
printTS("Solar", "DA")
printTS("Solar", "RT")
printTS("Wind", "DA")
printTS("Wind", "RT")

### MAIN LOOP ###
meanDailyWTrend <- NULL; meanDailySTrend <- NULL;
rep = 1; # dummy assignment
for (rep in 1:nRep) {
  # set the seed
  set.seed(rep)
  
  ### FIT THE MODEL ###
  source("runModel.R")
  
  # ~ mean scenario ~
  # compute the mean scenario (inc. lookahead)
  meanDailyWTrendtmp <- wModel$model$dailyTrend
  meanDailyWTrendtmp <- abind::abind(meanDailyWTrendtmp, wModel$model$dailyTrend[1:lookahead,], along = 1)

  meanDailySTrendtmp <- sModel$model$avgSky
  meanDailySTrendtmp <- abind::abind(meanDailySTrendtmp, sModel$model$avgSky[1:lookahead,], along = 1)

  # interpolate
  meanDailyWTrendtmp <- interpolate(simLength = simLength, lookahead = lookahead, wModel$numLoc, meanDailyWTrendtmp, 1)
  meanDailySTrendtmp <- interpolate(simLength = simLength, lookahead = lookahead, sModel$numLoc, meanDailySTrendtmp, 1)
  
  # clean up
  sunrise = (sModel$model$dayTime[1]-1-1)*4+2
  sunset = sModel$model$dayTime[ length(sModel$model$dayTime) ]*4
  meanDailySTrendtmp[-(sunrise:sunset),,] = 0
  
  meanDailyWTrendtmp <- meanDailyWTrendtmp[1:96,,1]
  meanDailySTrendtmp <- meanDailySTrendtmp[1:96,,1]
  
  # append
  meanDailyWTrend <- abind::abind(meanDailyWTrend, meanDailyWTrendtmp, along = 1)
  meanDailySTrend <- abind::abind(meanDailySTrend, meanDailySTrendtmp, along = 1)
  
  ### SIMULATE ###
  
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
  
  # continue with the next day
  today = today + 1
}

# append another mean-forecast (in case we look ahead)
meanDailyWTrend <- abind::abind(meanDailyWTrend, meanDailyWTrendtmp, along = 1)
meanDailySTrend <- abind::abind(meanDailySTrend, meanDailySTrendtmp, along = 1)

# print mean forecasts
write.csv(meanDailyWTrend, "Mean_Wind.csv")
write.csv(meanDailySTrend, "Mean_Solar.csv")

print("Scenario generation is completed.")

sink()
