# run.R
# Simulates the VAR-models

### Forecasts ###
wPaths = windSimulate(wModel = wModel, simLength = simLength, lookahead = lookahead, numScenarios = numScenarios)
sPaths = solarSimulate(sModel = sModel, simLength = simLength, lookahead = lookahead, numScenarios = numScenarios)

print("Simulation is completed")

### Interpolation for Obtaining Subhourly Data ###
wPaths <- interpolate(simLength = simLength, lookahead = lookahead, varmodel = wModel, varpaths = wPaths, numScenarios = numScenarios)
sPaths <- interpolate(simLength = simLength, lookahead = lookahead, varmodel = sModel, varpaths = sPaths, numScenarios = numScenarios)

# set supply to 0 during no-sun hours
sunrise = (sModel$model$dayTime[1])*4 - 2  # 1/2 is subtracted for a smoother sunrise
sunset = sModel$model$dayTime[ length(sModel$model$dayTime) ]*4
sPaths[-(sunrise:sunset),,] = 0

print("Interpolation is completed")