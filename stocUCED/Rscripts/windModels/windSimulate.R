windSimulate <- function(wModel, simLength = 24, lookahead = 3, numScenarios = 10) {
  # INPUT:
  #   - wModel    - wind statistical model
  #   - simLength - simulation length in hours
  #   - lookahead - # of lookahead periods in hours
  #   - numScenarios - number of scenarios to be simulated
  # OUTPUT:
  #   - A matrix with scenarios (time point x location x sample)
  
  # source('../tsUtilities/varSimulate.R')
  
  # Number of hours of data points to be computed
  nDays   <- simLength/24;
  nPointsPerDay <- 24*wModel$freq
  nPoints <- (simLength+lookahead)*wModel$freq
  
  # Simulate the residual time series
  simPaths <- varSimulate(varModel = wModel$model$est, sim.len  = nPoints, sim.freq = wModel$freq, 
                          numScenarios = numScenarios)
  
  # Normalize the values to be within [0,1] and then scale it using the generator capacity
  capacityMatrix <- replicate(n = numScenarios, t(replicate(n = nPoints, wModel$model$capacity)))
  simPaths <- simPaths/max(abs(simPaths))
  simPaths <- simPaths*capacityMatrix
  
  # extend the daily trends to cover the lookahead periods
  extDailyTrend <- matrix(nrow = (simLength+lookahead)*wModel$freq, ncol = wModel$numLoc)
  for (d in 1:nDays) {
    extDailyTrend[((d-1)*nPointsPerDay+1):(d*nPointsPerDay),] <- wModel$model$dailyTrend
  }
  if (lookahead > 24) {
    printf("lookahead must be at most 24 hours")
  } else if (lookahead > 0) {
    extDailyTrend[(nPointsPerDay*nDays+1):((simLength+lookahead)*wModel$freq),] <- wModel$model$dailyTrend[1:(lookahead*wModel$freq),]
  }
  
  # Loop across the days to align the clear sky model 
  samplePaths <- NULL
  for ( d in 1:nDays ) {
    tempTS <- replicate(n = numScenarios, expr = extDailyTrend[((d-1)*nPointsPerDay+1):(d*nPointsPerDay),]) + simPaths[((d-1)*nPointsPerDay+1):(d*nPointsPerDay),,]
    samplePaths <- abind::abind(samplePaths, tempTS, along = 1)
  }
  # lookahead..
  if (lookahead > 0) {
    tempTS <- replicate(n = numScenarios, expr = extDailyTrend[(nDays*nPointsPerDay+1):((simLength+lookahead)*wModel$freq),]) + simPaths[(nDays*nPointsPerDay+1):((simLength+lookahead)*wModel$freq),,]
    samplePaths <- abind::abind(samplePaths, tempTS, along = 1)
  }
  
  samplePaths[samplePaths < 0] <- 0;
  samplePaths[samplePaths > capacityMatrix] <- capacityMatrix[samplePaths > capacityMatrix]
  
  return(samplePaths)
}