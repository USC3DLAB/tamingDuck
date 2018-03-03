# run.R
# Fits the VAR-models

### Wind Model ###
tsData = prepareTSData(dataset = dataset, srcType = "Wind", prepType = "extract", forecastType = "DA", today = today, forecastWindow = forecastWindow)
wModel = windVARModel(tsData, lag.max = 5, infocrit = "SC")

### Solar Model ###
tsData = prepareTSData(dataset = dataset, srcType = "Solar", prepType = "extract", forecastType = "DA", today = today, forecastWindow = forecastWindow)
sModel = solarVARModel(tsData, lag.max = 5, infocrit = "SC", clearSky = TRUE)

print("Modeling is completed")