interpolate <- function(simLength, lookahead, varmodel, varpaths, numScenarios) {
  # Thx to: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
  library(parallel)
  nThreads <- detectCores()
  
  cl <- makeCluster(nThreads)
  
  subhourlyPaths <- array(0, dim = c((simLength+lookahead)*4, varmodel$numLoc, numScenarios) )
  for (gen in 1:varmodel$numLoc) {  # for every generator, interpolate the scenarios
    
    # this interpolation assumes the hourly data was recorded at HH:30.
    tmp <- parLapply(cl, 
                     (1:numScenarios), 
                     gen = gen, 
                     varpaths = varpaths, 
                     totLength = simLength+lookahead,
                     desiredFreq = 4,
                     function(j, gen, varpaths, totLength, desiredFreq) {
                       tmp <- spline(x = 1:totLength, y = varpaths[,gen,j], xout = seq(from = 1, to = (totLength+1-1e-6), by = 0.25))
                       tmp$y[ tmp$y<0 ] <- 0
                       return (tmp$y)
                     })
    
    # combine the results
    for (t in 1:((simLength+lookahead)*4/wModel$freq)) {
      for (s in 1:numScenarios) {
        subhourlyPaths[t,gen,s] = tmp[[s]][t]
      }
    }
  }
  stopCluster(cl)
  
  return(subhourlyPaths)
}