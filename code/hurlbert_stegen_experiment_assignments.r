# Assigning sims to scenarios and experiments within the Hurlbert & Stegen model
params = read.table('clipboard', header = T, sep = '\t', stringsAsFactors = F)

params$scenario[params$sim.id >= 6245 & params$energy.gradient=='on'] = 'K'
params$scenario[params$sim.id >= 6245 & params$specn.gradient == 'on'] = 'S'
params$scenario[params$sim.id >= 6245 & params$disturb_frequency > 0] = 'D'
params$scenario[params$sim.id >= 6245 & params$carry.cap == 'off'] = 'N'

params$env[params$sim.id >= 6245 & params$w == 7] = 'low'
params$env[params$sim.id >= 6245 & params$w == 5] = 'medium'
params$env[params$sim.id >= 6245 & params$w == 3 & params$sigma_E == 1 & params$beta == 1e-4 & params$alpha == 1e-6] = 'high'

params$nic[params$sim.id >= 6245 & params$sigma_E == 9] = 'low'
params$nic[params$sim.id >= 6245 & params$sigma_E == 5] = 'medium'
params$nic[params$sim.id >= 6245 & params$sigma_E == 1 & params$w == 3 & params$beta == 1e-4 & params$alpha == 1e-6] = 'high'

params$mut[params$sim.id >= 6245 & params$alpha == 1e-7] = 'low'
params$mut[params$sim.id >= 6245 & params$alpha == 1e-6 & params$sigma_E == 1 & params$beta == 1e-4 & params$w == 3] = 'medium'
params$mut[params$sim.id >= 6245 & params$alpha == 1e-5] = 'high'

params$dis[params$sim.id >= 6245 & params$beta == 1e-5] = 'low'
params$dis[params$sim.id >= 6245 & params$beta == 1e-4 & params$sigma_E == 1 & params$w == 3 & params$alpha == 1e-6] = 'medium'
params$dis[params$sim.id >= 6245 & params$beta == 1e-3] = 'high'


write.csv(params, 'parameters/hs_parameters.csv', row.names = F)
