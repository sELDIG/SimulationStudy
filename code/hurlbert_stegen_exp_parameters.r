# Generate simulation model parameter sets uniformly sampling relevant parameter space

paramLimits = data.frame(parameter = c('w', 'alpha', 'beta', 'sigma_E'),
                         min = c(1, -8, -7, 1),
                         max = c(11, -4, -2, 11))

sampleParams = function(paramLimits, n = 1000, seed = 1, startingSimID = 1) {
  set.seed(seed)
  
  output = data.frame(simID = startingSimID:(startingSimID + n - 1))
 
  for (i in 1:nrow(paramLimits)) {
    output = cbind(output, 
                   data.frame(runif(n)*(paramLimits$max[i] - paramLimits$min[i]) + paramLimits$min[i]))
  }
  names(output)[2:ncol(output)] = as.character(paramLimits$parameter)

  return(output)
}

## Generate SENC simulation matrix for each of 4 scenarios
randParams = sampleParams(paramLimits, n = 500, startingSimID = 6501)
randParams$alpha10 = 10^randParams$alpha
randParams$beta10 = 10^randParams$betaa

# 1) Energy gradient
EGParams = data.frame(sim.id = randParams$simID,
                      status = 'to.run',
                      reg.of.origin = 'tropical', 
                      w = randParams$w,
                      alpha = randParams$alpha10,
                      beta = randParams$beta10,
                      sigma_E = randParams$sigma_E,
                      carry.cap = 'on',
                      energy.gradient = 'on',
                      max.K = 10000,
                      num.of.bins = 11,
                      max.time = 30000,
                      max.richness = 5000,
                      replicate = 1,
                      disturb_frequency = 0,
                      temperate_disturb_intensity = 0,
                      tropical_disturb_intensity = 0,
                      specn.gradient = 'off',
                      specn.factor = NA,
                      gamma = 0.1)


# 2) Speciation gradient
SGParams = data.frame(sim.id = 7001:7500,
                      status = 'to.run',
                      reg.of.origin = 'tropical', 
                      w = randParams$w,
                      alpha = randParams$alpha10,
                      beta = randParams$beta10,
                      sigma_E = randParams$sigma_E,
                      carry.cap = 'on',
                      energy.gradient = 'off',
                      max.K = 10000,
                      num.of.bins = 11,
                      max.time = 30000,
                      max.richness = 5000,
                      replicate = 1,
                      disturb_frequency = 0,
                      temperate_disturb_intensity = 0,
                      tropical_disturb_intensity = 0,
                      specn.gradient = 'on',
                      specn.factor = 10,
                      gamma = 0.1)


# 3) Disturbance gradient
DGParams = data.frame(sim.id = 7501:8000,
                      status = 'to.run',
                      reg.of.origin = 'tropical', 
                      w = randParams$w,
                      alpha = randParams$alpha10,
                      beta = randParams$beta10,
                      sigma_E = randParams$sigma_E,
                      carry.cap = 'on',
                      energy.gradient = 'off',
                      max.K = 10000,
                      num.of.bins = 11,
                      max.time = 30000,
                      max.richness = 5000,
                      replicate = 1,
                      disturb_frequency = 100,
                      temperate_disturb_intensity = 0.99,
                      tropical_disturb_intensity = 0.75,
                      specn.gradient = 'off',
                      specn.factor = NA,
                      gamma = 0.1)


# 4) Time gradient
TGParams = data.frame(sim.id = 8001:8500,
                      status = 'to.run',
                      reg.of.origin = 'tropical', 
                      w = randParams$w,
                      alpha = randParams$alpha10,
                      beta = randParams$beta10,
                      sigma_E = randParams$sigma_E,
                      carry.cap = 'off',
                      energy.gradient = 'off',
                      max.K = 10000,
                      num.of.bins = 11,
                      max.time = 30000,
                      max.richness = 5000,
                      replicate = 1,
                      disturb_frequency = 0,
                      temperate_disturb_intensity = 0,
                      tropical_disturb_intensity = 0,
                      specn.gradient = 'off',
                      specn.factor = NA,
                      gamma = 0.1)

