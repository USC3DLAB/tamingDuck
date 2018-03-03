varSimulate <- function(varModel, sim.len, sim.freq, numScenarios) {
  # Input:
  #   varModel            = input of class attribute ‘varest’
  #   (sim.len, sim.freq) = simulation length (in hours) and frequency (number of data points per hour)
  #   numScenarios        = number of scenarios
  # Output:
  #   varSimulate         = a matrix with simulated outcome (each row is a new realization)
  
  # Thx to: http://gforge.se/2015/02/how-to-go-parallel-in-r-basics-tips/
  library(parallel)
  nThreads <- detectCores()
  
  simLength <- sim.len*sim.freq;
  # prepare to simulate scenarios from selected model.
  coeff <- NULL;
  for (i in 1:varModel$K) {
    coeff <- abind(coeff, varModel$varresult[[i]]$coefficients, along = 2);
  }
  
  # setting model order
  p <- varModel$p; q <- 0; skip <- 200; arlags <- 1:p; malags <- NULL;
  
  # setting the constant term
  cnst <- coeff[as.integer(varModel$K*p)+1,];
  
  # setting the AR lag-polynomial 
  Apoly <- t(coeff[1:as.integer(varModel$K*p),]);
  
  # setting the MA lag-polynomial 
  Bpoly <- NULL;
  
  # setting covariance to identify-matrix
  Spoly <- diag(x = 1, varModel$K);
  
  # simulate sample paths
  cl <- makeCluster(nThreads)
  samplePaths <- parLapply(cl, 1:numScenarios, 
                       nobs = simLength, arlags = arlags, malags = NULL, cnst = cnst,
                       phi = Apoly, theta = Bpoly, sigma = Spoly, skip = skip,
                       function(s, nobs, arlags, malags, cnst, phi, theta, sigma, skip) {
                         sample = MTS::VARMAsim(nobs = simLength, arlags = arlags, malags = malags, cnst = cnst, phi = Apoly, theta = Bpoly, sigma = Spoly, skip = skip)
                         return (sample$series)
                       })
  stopCluster(cl)
  samplePaths <- array(unlist(samplePaths), dim = c(nrow(samplePaths[[1]]), ncol(samplePaths[[1]]), length(samplePaths)))

  return(samplePaths)
}