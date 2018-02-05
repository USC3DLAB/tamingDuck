varModelFit <- function(dataFolder, fileName, dataFrequency, header, rowNames, infocrit) {
  require(vars); require(MTS); 
  
  # Load and preprocess time series 
  tsData <- preprocessNRELData(dataFolder, fileName , dataFrequency, header, rowNames)

  # Select model order according to chosen information criteria and estimate the coefficients of the chosen model.
  varsimest <- vars::VAR(tsData$detrend, type = "const", lag.max = 3, ic = infocrit);
  
  return(list(ts = tsData, model = varsimest));
}