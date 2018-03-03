prepareTSData <- function(dataset, srcType, prepType = "extract", forecastType = "DA",
                          today = 15, forecastWindow = 14,
                          decomposeBy = 'monthly', identifier = 'January', 
                          ts = NULL, cap = NULL) {
  # Input: 
  #   - dataset: the whole dataset
  #   - srcType: which source are we preparing for model-fitting? ("Solar", "Wind", ...)
  #   - prepType: "extract" extracts "forecastWindow" days before "today" from the assoc. time-series
  #               "decompose" extracts monthly, seasonally, daily, weekly data 
  #   - ts / cap: you can provide non-NULL arguments if you like to append
  # Output:
  #   - Extracted time-series in the form of a matrix where each row corresponds to one day's data. 
  
  # Read from data set corresponding to solar generators
  for (d in 1:length(dataset)) {
    if ( dataset[[d]]$srcType == srcType) {
      if (prepType == "extract") {
        ts <- abind::abind(ts, extractTSData(inputTS = dataset[[d]]$ts[[forecastType]], 
                                             freq = dataset[[d]]$freq[[forecastType]],
                                             today = today, forecastWindow = forecastWindow), along = 3);
      }
      else if (prepType == "decompose") {
        ts <- abind::abind(ts, decomposeTSData(inputTS = dataset[[d]]$ts[[forecastType]], 
                                               freq = dataset[[d]]$freq[[forecastType]], 
                                               decomposeBy = decomposeBy, identifier = identifier), along = 3)
      }
      else {
        print("Error: Unknown prepType in prepareTSData.R")
      }
      cap <- abind::abind(cap, dataset[[d]]$capacity)
      freq <- dataset[[d]]$freq[[forecastType]]
    }
  }
  ts <- aperm(ts, c(2,3,1))
  
  return(list(ts = ts, cap = cap, freq = freq))
}