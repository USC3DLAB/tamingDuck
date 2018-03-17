# run.R
# Simulates the VAR-models

### Forecasts ###
wPaths = windSimulate(wModel = wModel, simLength = simLength, lookahead = lookahead, numScenarios = numScenarios)
sPaths = solarSimulate(sModel = sModel, simLength = simLength, lookahead = lookahead, numScenarios = numScenarios, clearSky = FALSE)

print("Simulation is completed")

### Interpolation for Obtaining Subhourly Data ###
wPaths <- interpolate(simLength = simLength, lookahead = lookahead, nComponents = wModel$numLoc, varpaths = wPaths, numScenarios = numScenarios)
sPaths <- interpolate(simLength = simLength, lookahead = lookahead, nComponents = sModel$numLoc, varpaths = sPaths, numScenarios = numScenarios)

# set supply to 0 during no-sun hours
sunrise = (sModel$model$dayTime[1]-1-1)*4+2
sunset = sModel$model$dayTime[ length(sModel$model$dayTime) ]*4
sPaths[-(sunrise:sunset),,] = 0

print("Interpolation is completed")