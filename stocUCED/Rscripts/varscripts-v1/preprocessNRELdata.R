preprocessNRELData <- function(dataFolder, fileName, dataFrequency, header, rowNames) {
  # Read the time series data and preprocess it
  # Input: 
  #   dataFolder    = "/home/gjharsha/Documents/workspace/tamingDuck/stocUCED/datasets/3d-nrel118/Wind"
  #   filename      = "DA.csv";
  #   dataFrequency = Number of data points in one hour (1 or 4)

  require("abind");
  
  # Setup the full path for the files and read them.
  fullPath = sprintf('%s/%s', dataFolder, fileName);
  temp <- read.csv(file = fullPath, header = header, row.names = rowNames);
  
  # Data summary: number of locations and time periods (in hours)
  tsData = list(N = dim(temp)[2], T = dim(temp)[1]/dataFrequency, freq = dataFrequency)

  # Preprocess the data
  tsData$data <- temp; tsData$detrend <- NULL;
  tsData$dailyData <- NULL; tsData$dailyMean <- NULL; 
  for (l in 1:tsData$N) {
    # decompose the time series into daily time series
    tsData$dailyData <- abind(tsData$dailyData, t(matrix(data = tsData$data[,l], nrow = 24)), along = 3)
    
    # Compute the daily mean
    tsData$dailyMean <- abind(tsData$dailyMean, apply(tsData$dailyData[,,l], MARGIN = 2, FUN = mean), along = 2);
    
    M <- NULL
    # De-trend the data by removing the daily mean
    for ( d in 1:dim(tsData$dailyData)[1]) {
      M <- abind(M, tsData$dailyData[d,,l] - tsData$dailyMean[l], along = 2)
    }
    
    tsData$detrend <- abind(tsData$detrend, matrix(data = M, ncol = 1), along = 2); 
  }
  
  return(tsData);
}