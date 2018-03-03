extractTSData <- function(inputTS, freq, today = 15, forecastWindow = 14) {
  # Input: 
  #   - inputTS = A one-dimensional time series (assumed to be annual data)
  #   - freq = Frequency of data
  #   - today = current date index (days)
  #   - forecastWindow = forecasting will be done based on this many days
  # Output:
  #   - Extracted time-series in the form of a matrix where each row corresponds to one day's data. 
  
  today = (today-1) * freq * 24 + 1;
  forecastWindow = forecastWindow * freq * 24;
  
  # extract the selected indices from the input time-series
  selectedIndices = (today-forecastWindow):(today-1)
  inputTS <- inputTS[selectedIndices]
  
  # make the time-series daily
  tsData <- t(matrix(data = inputTS, nrow = freq*24))
  
  return(tsData)
}