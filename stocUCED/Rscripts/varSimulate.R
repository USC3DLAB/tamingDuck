varSimulate <- function(varModel, simLength, numScenarios, simFrequency, type) {
  
  simLength <- simLength*simFrequency;
  trend <- NULL
  for ( l in 1:varModel$ts$N ) {
    temp <- spline(x = 1:dim(varModel$ts$dailyMean)[1], y = varModel$ts$dailyMean[,l], n = simLength)
    temp$y[temp$y < 0] <- 0
    trend <- abind(trend, temp$y, along = 2)
  }
  
  samplePaths <- NULL;
  for (s in 1:numScenarios) {  
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
    windowSim <- MTS::VARMAsim(nobs = simLength, arlags = arlags, malags = NULL, cnst = cnst, phi = Apoly, theta = Bpoly, sigma = Spoly, skip = skip);
    
    oneSample <- windowSim$noises + trend;
    oneSample[oneSample < 0] <- 0
    if ( type == "Solar" )
      oneSample[trend == 0] <- 0;
    
    samplePaths <- abind(samplePaths, oneSample, along = 3);
  }
  
  return(samplePaths)
}