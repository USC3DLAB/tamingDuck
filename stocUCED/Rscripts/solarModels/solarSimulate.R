solarSimulate <- function(sModel, simLength = 24, lookahead = 3, numScenarios = 10) {
  # INPUT:
  #   - sModel    - solar statistical model
  #   - simLength - simulation length in hours
  #   - lookahead - # of lookahead periods in hours
  #   - numScenarios - number of scenarios to be simulated
  # OUTPUT:
  #   - A matrix with scenarios (time point x location x sample)
  
  # source('../tsUtilities/varSimulate.R')
  
  # Number of hours of data points to be computed
  nDays   <- simLength/24;
  nPoints <- nDays*length(sModel$model$dayTime)
  
  # Simulate the residual time series
  simPaths <- varSimulate(varModel = sModel$model$est, sim.len  = nPoints, sim.freq = sModel$freq, 
                             numScenarios = numScenarios)
  
  # Normalize the values to be within [0,1]
  if ( min(simPaths) < 0 )
    simPaths <- simPaths - min(simPaths)
  simPaths <- simPaths/max(simPaths)
  
  # Loop across the days to align the clear sky model 
  nDayHours <- length(sModel$model$dayTime)
  samplePaths <- NULL
  for ( d in 1:nDays ) {
    tempTS <- replicate(n = numScenarios, expr = sModel$model$clearSky);
    tempTS[sModel$model$dayTime,,] <- tempTS[sModel$model$dayTime,,]*simPaths[((d-1)*nDayHours+1):(d*nDayHours),,]
    samplePaths <- abind::abind(samplePaths, tempTS, along = 1)
  }
  
  if (lookahead >= sModel$model$dayTime[1]) {
    print("Error: I cannot lookahead another day-time, during solar forecasting")
  } else if (lookahead > 0) {
    for (l in 1:lookahead) {
      samplePaths <- abind::abind(samplePaths, matrix(data = 0, nrow = sModel$numLoc, ncol = numScenarios), along = 1)
    }
  }
  return(samplePaths)
}