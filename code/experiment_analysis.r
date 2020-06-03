# Join parameter files indicating high, medium, and low treatment levels to tree output.

library(dplyr)
library(nlme)

source('code/pcaFunctions.r')
source('code/experiment_analysis_functions.r')


experiments = data.frame(experiment = c('env', 'nic', 'dis', 'mut', 'tim'), 
                         phrase = c('environmental filtering', 'niche conservatism', 'disperal', 'mutation/speciation rate', 'time'))


# Read in output and join to parameter files
metrics = read.table('experiment_output.txt', header = T, sep = '\t', stringsAsFactors = FALSE)

models = unique(metrics$model)

joinedOutput = data.frame(model = NA, model2 = NA, simID = NA, S = NA, log10S = NA, tree.length = NA, PD = NA, Gamma = NA, 
                          Beta = NA, Colless = NA, Sackin = NA, Yule.PDA.ratio = NA, MRD = NA, 
                          VRD = NA, PSV = NA, mean.Iprime = NA, MPD = NA, VPD = NA, 
                          MGL_principal_eigenvalue = NA, MGL_asymmetry = NA,
                          MGL_peakedness = NA, MGL_eigengap = NA, nLTT_stat = NA,
                          env = NA, nic = NA, dis = NA, mut = NA, tim = NA)

for (m in models) {
  param = read.csv(paste('parameters/', m, '_parameters.csv', sep = ''), header = T, stringsAsFactors = FALSE)
  
  param$model2 = param$model
  param$model2[param$model == 'hs'] = paste(param$model, param$scenario, sep = '.')
  
  treatments = param[, names(param) %in% c('model', 'model2', 'simID', 'env', 'nic', 'dis', 'mut', 'tim')]
  
  # Add experiment columns if missing
  if(! 'env' %in% names(treatments)) treatments$env = NA
  if(! 'nic' %in% names(treatments)) treatments$nic = NA
  if(! 'dis' %in% names(treatments)) treatments$dis = NA
  if(! 'mut' %in% names(treatments)) treatments$mut = NA
  if(! 'tim' %in% names(treatments)) treatments$tim = NA
  
  # Order columns for consistency
  treatmentsOrdered = treatments %>% dplyr::select(model, model2, simID, env, nic, dis, mut, tim)
  
  # Join treatment columns to metrics output for a given model
  metricsSubset = filter(metrics, model == m) %>%
    left_join(treatmentsOrdered, by = c('model', 'simID')) %>%
    dplyr::select(model, model2, simID, S:tim)
  
  joinedOutput = rbind(joinedOutput, metricsSubset)
  
}

modelColors = data.frame(model2 = unique(joinedOutput$model2), color = colorSelection(length(unique(joinedOutput$model2))))
modelColors$model2 = as.character(modelColors$model2)
modelColors$color = as.character(modelColors$color)


expOutput = joinedOutput[-1, ] %>%
  mutate(envStandardized = toupper(substr(env, 1, 1)),
         nicStandardized = toupper(substr(nic, 1, 1)),
         disStandardized = toupper(substr(dis, 1, 1)),
         mutStandardized = toupper(substr(mut, 1, 1)),
         timStandardized = toupper(substr(tim, 1, 1)),
         envLevel = case_when(
           envStandardized == 'L' ~ 1,
           envStandardized == 'M' ~ 2,
           envStandardized == 'H' ~ 3),
         nicLevel = case_when(
           nicStandardized == 'L' ~ 1,
           nicStandardized == 'M' ~ 2,
           nicStandardized == 'H' ~ 3),
         disLevel = case_when(
           disStandardized == 'L' ~ 1,
           disStandardized == 'M' ~ 2,
           disStandardized == 'H' ~ 3),
         mutLevel = case_when(
           mutStandardized == 'L' ~ 1,
           mutStandardized == 'M' ~ 2,
           mutStandardized == 'H' ~ 3),
         timLevel = case_when(
           timStandardized == 'L' ~ 1,
           timStandardized == 'M' ~ 2,
           timStandardized == 'H' ~ 3)) %>%
  left_join(modelColors, by = 'model2')




## Plotting results

plotExperimentResults('env', expOutput)
plotExperimentResults('nic', expOutput)
plotExperimentResults('dis', expOutput)
plotExperimentResults('mut', expOutput)
plotExperimentResults('tim', expOutput)



# Correlation output
corrOutput = data.frame(model2 = NA, 
                    experiment = NA,
                    r.log10S = NA,
                    r.PD = NA,
                    r.Gamma = NA,
                    r.Beta = NA,
                    r.Colless = NA,
                    r.Sackin = NA,
                    r.Yule.PDA.ratio = NA,
                    r.MRD = NA,
                    r.VRD = NA,
                    r.PSV = NA,
                    r.mean.Iprime = NA,
                    r.MPD = NA,
                    r.VPD = NA,
                    r.nLTT_stat = NA)

for (m in unique(expOutput$model2)) {
  
  for (e in c('env', 'nic', 'dis', 'mut', 'tim')) {
    corrs = corrCalc(experiment = e, expOutput, mod = m)
  
    if (!is.na(corrs)) {
      corrOutput = rbind(corrOutput, corrs)  
    }
  }
}
corrOutput = corrOutput[-1, ] %>%
  left_join(modelColors, by = 'model2') %>%
  arrange(experiment, model2)


# Plotting histograms of correlations
pdf('figures/corr_histograms.pdf', height = 8, width = 10)
par(mar = c(4, 4, 0, 0), oma = c(0, 0, 3, 0), mgp = c(2.5, 1, 0), mfrow = c(4, 4))

for (experiment in c('env', 'nic', 'dis', 'mut', 'tim')) {
  for (met in names(corrOutput[3:ncol(corrOutput)])) {
      hist(corrOutput[corrOutput$experiment == experiment, met], xlab = met, yaxt = 'n',
           main = "", xlim = c(-1, 1), breaks = seq(-1, 1, by = .2), col = 'gray50')
      axis(2, at = 0:3, las = 1)
      abline(v = 0, col = 'red', lwd = 3)
  }
  mtext(paste(experiments$phrase[experiments$experiment == experiment], "experiment"), 3, outer = T, cex = 3)
  par(mar = c(4, 4, 0, 0), oma = c(0, 0, 3, 0), mgp = c(2.5, 1, 0), mfrow = c(4, 4))
}
dev.off()


# Plotting correlation coefficients by model
pdf('figures/corr_plots.pdf', height = 8, width = 10)
par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), mfrow = c(4, 4), cex.axis = 1.3)

for (exp in c('env', 'nic', 'dis', 'mut', 'tim')) {
  
  for (met in names(corrOutput[3:(ncol(corrOutput)-1)])) {
    tmp = filter(corrOutput, experiment == exp) %>% 
      arrange(get(met))
    plot(tmp[, met], 1:nrow(tmp), pch = 16, col = tmp$color, 
         cex = 2, xlim = c(-1, 1), xlab = '', yaxt = 'n', ylab = '')
    text(-1, nrow(tmp), met, cex = 1.5, adj = c(0, 1))
    abline(v = 0, col = 'black', lwd = 2)
  }
  # legend panel
  plot(1, 1, type = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', bty = 'n')
  points(rep(0.8, nrow(tmp)), seq(0.6, 1.4, length.out = nrow(tmp)),
         pch = 16, col = tmp$color, cex = 2)
  text(rep(1, nrow(tmp)), seq(0.6, 1.4, length.out = nrow(tmp)), tmp$model2, cex = 1.5)
  
  
  mtext(paste(experiments$phrase[experiments$experiment == exp], "experiment"), 3, outer = T, cex = 2, line = 1)
  mtext("correlation coefficient", 1, outer = T, cex = 1.5, line = 0, at = .3)
  
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), mfrow = c(4, 4))
}
dev.off()




# Mixed effects model output
metrics = names(expOutput)[c(5, 7:18, 23)]

MEoutput = data.frame(experiment = NA, 
                      metric = NA,
                      model = NA,
                      slope = NA,
                      intercept = NA,
                      ngroups = NA,
                      totalObs = NA,
                      global.slope.sd = NA,
                      global.int.sd = NA)

for (experiment in c('env', 'nic', 'dis', 'mut', 'tim')) {
  
  # Check whether there are multiple models available to run lmer for a given experiment
  if(length(unique(expOutput$model2[!is.na(expOutput[[paste(experiment, "Level", sep = "")]])])) >= 2) {
    
    for (met in metrics) {
      
      lme.mod = lmer(expOutput[[met]] ~ expOutput[[paste(experiment, "Level", sep = "")]] + 
                       (1 + expOutput[[paste(experiment, "Level", sep = "")]] | model2), 
                     data = expOutput, REML = FALSE)
      
      lmeOut = extractLMEoutput(lme.mod, experiment, met)
      
      MEoutput = rbind(MEoutput, lmeOut)
    }
    
  }
}
