solarVARModel <- function(tsData, freq, lag.max = 3, infocrit = 'SC', clearSky = TRUE) {
  # Setup and fit a VAR model for the solar dataSet
  #   - tsData = time-series related data
  #   - freq = dataset frequency
  #   - lag.max     = maximum lag order for model estimation, 
  #   - infocrit    = information criteria used for model order selection
  
  ts = tsData$ts; cap = tsData$cap; freq = tsData$freq;
  
  # Data parameters
  numDataPoints <- dim(ts)[1]; numLoc <- dim(ts)[2]; numDays <- dim(ts)[3]
  
  # Define the clear sky generation at a location to be the maximum generation across all days and locations
  dayTime <- which(x = apply(X = ts, MARGIN = 1, FUN = max) > 0)
  
  # Normalize the using the clear sky model 
  clearSkyGeneration <- NULL; normalizedTS <- NULL; avgSkyGeneration <- NULL;
  for ( n in 1:numLoc ) {
    # Find the clear-sky generation for the location
    clearSkyGeneration <- abind::abind(clearSkyGeneration, apply(X = ts[,n,], MARGIN = 1, FUN = max), along = 2);
    
    # Find the average-sky generation for the location
    avgSkyGeneration <- abind::abind(avgSkyGeneration, apply(X = ts[,n,], MARGIN = 1, FUN = mean), along = 2);
    
    if (clearSky) {
      # replicate the clear sky generation to apply for all days in the horizon
      scr <- replicate(n = numDays, clearSkyGeneration[,n])
      
      # Normalize the data using the clear sky generation (replace NaN with zeros)
      tempTS <- ts[,n,]/scr;
      tempTS[is.nan(tempTS)] <- 0;
      
    } 
    else {
      # replicate the avg sky generation to apply for all days in the horizon
      scr <- replicate(n = numDays, avgSkyGeneration[,n])
      
      # Normalize the data using the avg sky generation
      tempTS <- ts[,n,] - avgSkyGeneration[,n]
      tempTS <- tempTS/cap[n]
    }
    
    # Save the results only for the day time
    tempTS <- tempTS[dayTime,]
    
    # Resize the data as a single normalized time series
    normalizedTS <- abind::abind(normalizedTS, as.vector(tempTS), along = 2)
  }
  
  # Use the normalized data to fit a VAR model
  varsimest <- vars::VAR(normalizedTS, type = "const", lag.max = lag.max, ic = infocrit);
  
  # check for NA coefficients
  for (name in names(varsimest$varresult)) {
    if (anyNA(varsimest$varresult[[name]]$coefficients)) {
      print("Error! The VAR model cannot estimate regression coefficients!")
    }
  }
  
  return(list(ts = ts, numLoc = numLoc, freq = freq,
              model = list(clearSky = clearSkyGeneration, avgSky = avgSkyGeneration, capacity = cap, dayTime = dayTime, est = varsimest)))
}