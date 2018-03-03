interpolate <- function(simLength, lookahead, nComponents, varpaths, numScenarios) {
  # Thx to: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
  library(parallel)
  nThreads <- detectCores()
  
  # if a matrix, instead of 3D array is provided, do the conversion
  if (numScenarios == 1 && length(dim(varpaths)) == 2) { 
    varpaths = array(varpaths, dim = c(dim(varpaths), 1))
  }
  
  cl <- makeCluster(nThreads)
  
  subhourlyPaths <- NULL
  for (gen in 1:nComponents) {  # for every generator, interpolate the scenarios
    
    # this interpolation assumes the hourly data was recorded at HH:00.
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
    subhourlyPaths <- abind::abind(subhourlyPaths, array(unlist(tmp), c(length(tmp[[1]]), length(tmp))), along = 3)
#    for (t in 1:((simLength+lookahead)*4)) {
#      for (s in 1:numScenarios) {
#        subhourlyPaths[t,gen,s] = tmp[[s]][t]
#      }
#    }
  }
  stopCluster(cl)
  
  subhourlyPaths <- aperm(subhourlyPaths, c(1,3,2))
  
  return(subhourlyPaths)
}