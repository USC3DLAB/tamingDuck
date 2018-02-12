# print scenarios
for (i in 1:numScenarios) {
  write.csv(scenarios[,,i], sprintf("sim_rep%d_scen%d.csv", rep+1, i))
}
