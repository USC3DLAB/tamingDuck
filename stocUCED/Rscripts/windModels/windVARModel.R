windVARModel <- function(dataSet, decomposeBy = 'monthly', identifier = 'February', forecastType = 'DA',
                         lag.max = 3, infocrit = 'SC') {
  # Setup and fit a VAR model for the solar dataSet
  #   - dataSet     = standard dataSet
  #   - decomposeBy = 'seasonal', 'monthly', 'weekly', or 'daily'
  #   - identifier  = ('summer', 'autumn', 'winter', 'spring') or name of the month
  #   - forecastType = ('DA', 'RT')
  #   - lag.max     = maximum lag order for model estimation, 
  #   - infocrit    = information criteria used for model order selection)
  
  # source('../tsUtilities/decomposeTSData.R')
  
  # Read from data set corresponding to wind generators
  ts <- NULL; cap <- NULL
  for (d in 1:length(dataSet)) {
    if ( dataSet[[d]]$srcType == 'Wind' ) {
      ts <- abind(ts, decomposeTSData(inputTS = dataSet[[d]]$ts[[forecastType]], 
                                      freq = dataSet[[d]]$freq[[forecastType]], 
                                      decomposeBy = decomposeBy, identifier = identifier), along = 3)
      cap <- c(cap, dataSet[[d]]$capacity)
      freq <- dataSet[[d]]$freq[[forecastType]]
    }
  }
  ts <- aperm(ts, c(2,3,1))
  
  # Data parameters
  numDataPoints <- dim(ts)[1]; numLoc <- dim(ts)[2]; numDays <- dim(ts)[3]
  
  # Normalize the data by removing the daily trend 
  dailyTrend <- NULL; normalizedTS <- NULL
  for ( n in 1:numLoc ) {
    # Compute trend as the mean of the nornmalized time series across all the days
    # We normalize the time series by the capacity of the generator
    dailyTrend <- abind::abind(dailyTrend, apply(X = ts[,n,], MARGIN = 1, FUN = mean), along = 2)
    
    # Remove the trend
    tempTS <- ts[,n,] - replicate(n = numDays, expr = dailyTrend[,n])
    
    # normalize the dataset
    tempTS <- tempTS/cap[n]
    # tempTS <- (tempTS - mean(as.vector(tempTS)))/var(as.vector(tempTS))
    
    # Resize the data as the normalize time series
    normalizedTS <- abind::abind(normalizedTS, as.vector(tempTS), along = 2)
  }
  
  # Use the normalized data to fit a VAR model
  varsimest <- vars::VAR(normalizedTS, type = "const", lag.max = lag.max, ic = infocrit);
  
  return(list(ts = ts, numLoc = numLoc, freq = freq, model = list(dailyTrend = dailyTrend, capacity = cap, est = varsimest)))
}