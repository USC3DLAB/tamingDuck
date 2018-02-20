varParallelSimulate <- function(varModel, simLength, lookaheadPeriods, simFrequency, numScenarios, type) {

  # Thx to: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
  library(parallel)
  
  # Calculate the number of threads
  nbThreads <- detectCores() - 1
  
  source("varParallelSimulateSub.R")
  
  simLength <- simLength*simFrequency;
  trend <- matrix(nrow = (simLength+lookaheadPeriods), ncol = varModel$ts$N);
  for ( l in 1:varModel$ts$N ) {
    temp <- spline(x = 1:dim(varModel$ts$dailyMean)[1], y = varModel$ts$dailyMean[,l], n = simLength)
    temp$y[temp$y < 0] <- 0
    trend[1:simLength, l] <- temp$y;
    trend[(simLength+1):(simLength+lookaheadPeriods), l] <- temp$y[1:(lookaheadPeriods)];
  }
  
  # Initiate cluster
  cl <- makeCluster(nbThreads)
  
  # Execute simulations
  hop <- parLapply(cl, 1:numScenarios, varParallelSimulateSub,
                varModel = varModel, 
                simLength = simLength,
                lookaheadPeriods = lookaheadPeriods,
                type = type,
                trend = trend)
  
  # Stop the cluster
  stopCluster(cl)
  
  # merge the scenarios
  samplePaths <- NULL;
  for (s in 1:numScenarios) {
    samplePaths <- abind(samplePaths, hop[[s]], along = 3);
  }
  
  return(samplePaths)
}
