varParallelSimulateSub <- function(s, varModel, simLength, lookaheadPeriods, type, trend) {
  library("abind");
  
  # set seed 
  # print(s);
  # set.seed(s);
  
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
  
  # done
  return(oneSample)
}
