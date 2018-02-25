windSimulate <- function(wModel, simLength = 24, numScenarios = 10) {
  # INPUT:
  #   - wModel    - wind statistical model
  #   - simLength - simulation length in hours
  #   - numScenarios - number of scenarios to be simulated
  # OUTPUT:
  #   - A matrix with scenarios (time point x location x sample)
  
  source('../tsUtilities/varSimulate.R')
  
  # Number of hours of data points to be computed
  nDays   <- simLength/24;
  nPointsPerDay <- 24*wModel$freq
  nPoints <- nDays*nPointsPerDay
  
  # Simulate the residual time series
  simPaths <- varSimulate(varModel = wModel$model$est, sim.len  = nPoints, sim.freq = wModel$freq, 
                          numScenarios = numScenarios)
  
  # Normalize the values to be within [0,1] and then scale it using the generator capacity
  capacityMatrix <- replicate(n = numScenarios, t(replicate(n = nPoints, wModel$model$capacity)))
  simPaths <- simPaths/max(abs(simPaths))
  simPaths <- simPaths*capacityMatrix
  
  # Loop across the days to align the clear sky model 
  samplePaths <- NULL
  for ( d in 1:nDays ) {
    tempTS <- replicate(n = numScenarios, expr = wModel$model$dailyTrend) + simPaths[((d-1)*nPointsPerDay+1):(d*nPointsPerDay),,]
    samplePaths <- abind::abind(samplePaths, tempTS, along = 2)
  }
  
  samplePaths[samplePaths < 0] <- 0;
  samplePaths[samplePaths > capacityMatrix] <- capacityMatrix[samplePaths > capacityMatrix]
  
  return(samplePaths)
}