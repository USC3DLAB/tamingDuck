varSimulate <- function(varModel, simLength, lookaheadPeriods, simFrequency, numScenarios, type) {
  
  simLength <- simLength*simFrequency;
  trend <- matrix(nrow = (simLength+lookaheadPeriods), ncol = varModel$ts$N);
  for ( l in 1:varModel$ts$N ) {
    temp <- spline(x = 1:dim(varModel$ts$dailyMean)[1], y = varModel$ts$dailyMean[,l], n = simLength)
    temp$y[temp$y < 0] <- 0
    trend[1:simLength, l] <- temp$y;
    trend[(simLength+1):(simLength+lookaheadPeriods), l] <- temp$y[1:(lookaheadPeriods)];
  }
  
  samplePaths <- NULL;
  for (s in 1:numScenarios) {
	  
	  if (s %% 100 == 0) {
		  print(s);
	  }
	
	  residualPath <- NULL;
    
    # prepare to simulate scenarios from selected model.
    coeff <- NULL;
    for (i in 1:varModel$ts$N) {
      coeff <- abind(coeff, varModel$model$varresult[[i]]$coefficients, along = 2);
    }
    
    # setting model order
    p <- varModel$model$p; q <- 0; skip <- 200; 
    arlags <- 1:p; malags <- NULL;
    
    # setting the constant term
    cnst <- coeff[as.integer(varModel$ts$N*p)+1,];
    
    # setting the AR lag-polynomial 
    Apoly <- t(coeff[1:as.integer(varModel$ts$N*p),]);
    
    # setting the MA lag-polynomial 
    Bpoly <- NULL;
    
    # setting covariance to identify-matrix
    Spoly <- diag(x = 1, varModel$ts$N);
    
    # simulate scenarios using VARMAsim
    windowSim <- MTS::VARMAsim(nobs = (simLength+lookaheadPeriods), arlags = arlags, malags = NULL, cnst = cnst, phi = Apoly, theta = Bpoly, sigma = Spoly, skip = skip);
    
    oneSample <- windowSim$noises + trend;
    oneSample[oneSample < 0] <- 0
	  if ( type == "Solar" ) {
		  for (g in 1:varModel$ts$N) {
		    oneSample[ trend[,g]/max(trend[,g]) < 1e-2 , g] = 0; 
		  }
	  }
	
    samplePaths <- abind(samplePaths, oneSample, along = 3);
  }
  
  return(samplePaths)
}
