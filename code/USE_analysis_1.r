
# Load libraries
library(gsheet)
library(dplyr)
library(stringr)
library(ape)
library(geiger)
library(lessR)

source('code/pcaFunctions.r')
source('code/experiment_analysis_functions.r')
source('code/treeMetrics.r')

# Read in dataframe of tree metrics and simulation parameters
processDFmetrics = read.csv('experiments/uniform_sampling_experiment/process_parameter_values_and_tree_metrics.csv')

paramKey <- read.csv('experiments/uniform_sampling_experiment/simulation_parameters_key.csv')








#########################################################################################
# Generate correlations between all parameters relevant to experiments and tree metrics

experiments = data.frame(experiment = c('env', 'nic', 'dis', 'mut', 'com'), 
                         phrase = c('environmental filtering', 'niche conservatism', 'dispersal', 
                                    'mutation/speciation rate', 'competition'))


corrOutput = data.frame(model = character(),
                        experiment = character(),
                        parameter = character(),
                        metric = character(),
                        r = double(),
                        r.L95 = double(),
                        r.U95 = double())

metriccolumns = processDFmetrics %>%
  dplyr::select(S:VPD, nLTT_stat:PC6) %>% names()

for (e in experiments$experiment) {
  
  relevantModels = filter(paramKey, experiment == e) %>% 
    distinct(model) %>% unlist()
  
  for (mod in relevantModels) {
    
    # Check for new model-experiment combinations to run
    if (!paste(mod, e) %in% paste(corrOutput$model, corrOutput$experiment)) {
      
      corDF = corrCalcUSE(processDFmetrics, e, mod, 
                          metrics = metriccolumns)
      
      if (class(corDF) == 'data.frame') {
        corrOutput = rbind(corrOutput, corDF)
      }
      
    }
    
  }
  
}


models = unique(paramKey$model)  #### Change this back to paramKey once gen is fixed
modelColors = data.frame(model = models, color = colorSelection(length(models)))
modelColors$model = as.character(modelColors$model)
modelColors$color = as.character(modelColors$color)

corrOutput2 = left_join(corrOutput, modelColors, by = 'model')

write.csv(corrOutput2, "experiments/uniform_sampling_experiment/treeMetric_parameter_correlation_output.csv", row.names = F)




# Summarizing best correlations for each model-experiment

maxcorr = corrOutput2 %>% 
  filter(metric != "S", !grepl("PC", metric)) %>% 
  group_by(model, experiment) %>% 
  summarize(metric2 = metric[abs(r) == max(abs(r))], 
            param = parameter[abs(r) == max(abs(r))],
            r = r[abs(r) == max(abs(r))], 
            r.L95 = r.L95[metric == metric2 & parameter == param], 
            r.U95 = r.U95[metric == metric2 & parameter == param]) %>%
  rename(metric = metric2,
         parameter = param) %>%
  left_join(data.frame(experiment = experiments$experiment, 
                       color = colorSelection(length(experiments$experiment)), 
                       ypos = 1:length(experiments$experiment)), by = 'experiment')




####################################
# Plotting the maximum parameter-metric correlation for each experiment in each model

pdf(paste0('figures/USE_max_corrs_by_model_and_experiment_', Sys.Date(), '.pdf'), height = 8, width = 10)
par(mar = c(2, 5, 1, 1), mfrow = c(4, 4))
for (m in unique(maxcorr$model)) {
  
  tmp = filter(maxcorr, model == m)
  
  plot(tmp$r, tmp$ypos, pch = 18, col = tmp$color, 
       cex = 1.8, xlim = c(-1.3, 1.3), xlab = '', yaxt = 'n', ylab = '', ylim = c(0.5, 6.5), cex.lab = .8)
  abline(v = 0, col = 'black', lwd = 1.2)
  segments(tmp$r.L95, tmp$ypos, tmp$r.U95, tmp$ypos, col = tmp$color, lwd = 2.5)
  text(-1.3, 6.5, m, cex = 2, adj = c(0, 1))
  
  text(tmp$r, tmp$ypos + .3, paste(tmp$parameter, "-", tmp$metric), cex = 0.75)
  
  mtext(tmp$experiment, 2, at = tmp$ypos, las = 1, line = 1, cex = 1, col = tmp$color)
}
dev.off()








# Plotting correlation coefficients by model for individual tree metrics
pdf(paste0('figures/USE_corr_plots_', Sys.Date(), '.pdf'), height = 8, width = 10)

for (exp in c('env', 'nic', 'dis', 'mut', 'com')) {
  
  layout(matrix(c(1:15, 15), nrow = 4, byrow = T))
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), cex.axis = 1.3)
  
  for (met in c('log10S', 'PD', 'Gamma', 'Beta', 'Colless', 'Sackin', 'Yule.PDA.ratio', 
                'MRD', 'VRD', 'PSV', 'mean.Iprime', 'MPD', 'VPD', 'nLTT_stat')) {
    
    tmp = filter(corrOutput2, experiment == exp, metric == met)
    
    plot(tmp$r, 1:nrow(tmp), pch = 18, col = tmp$color, 
         cex = 1.8, xlim = c(-1, 1), xlab = '', yaxt = 'n', ylab = '', ylim = c(0.5, 1.1*nrow(tmp)))
    abline(v = 0, col = 'black', lwd = 2)
    segments(tmp$r.L95, 1:nrow(tmp), tmp$r.U95, 1:nrow(tmp), col = tmp$color, lwd = 2.5)
    text(-1, nrow(tmp), met, cex = 1.5, adj = c(0, 1))
    
  }
  # legend panel
  plot(1, 1, type = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', bty = 'n')
  
  numMods = nrow(tmp)
  if (numMods <= 6) {
    points(rep(0.58, nrow(tmp)), seq(0.62, 1.38, length.out = nrow(tmp)),
           pch = 18, col = tmp$color, cex = 2)
    text(rep(.61, nrow(tmp)), seq(0.62, 1.38, length.out = nrow(tmp)), 
         paste(tmp$model, '-', tmp$parameter), cex = 1.3, adj = 0)
    
  } else {
    firstHalf = round(numMods/2)
    
    points(rep(0.58, firstHalf), seq(0.62, 1.38, length.out = firstHalf),
           pch = 18, col = tmp$color[1:firstHalf], cex = 2)
    text(rep(.61, firstHalf), seq(0.62, 1.38, length.out = firstHalf), 
         paste(tmp$model[1:firstHalf], '-', tmp$parameter[1:firstHalf]), cex = 1.1, adj = 0)
    
    points(rep(1.08, (numMods - firstHalf)), seq(0.62, 1.38, length.out = (numMods - firstHalf)),
           pch = 18, col = tmp$color[(firstHalf + 1):numMods], cex = 2)
    text(rep(1.12, (numMods - firstHalf)), seq(0.62, 1.38, length.out = (numMods - firstHalf)), 
         paste(tmp$model[(firstHalf + 1):numMods], '-', tmp$parameter[(firstHalf + 1):numMods]), cex = 1.1, adj = 0)
    
  }
  
  mtext(paste(experiments$phrase[experiments$experiment == exp], "experiment"), 3, outer = T, cex = 2, line = 1)
  mtext("correlation coefficient", 1, outer = T, cex = 1.5, line = 0, at = .3)
  
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), mfrow = c(4, 4))
}
dev.off()



# Plotting correlation coefficients by model for tree metric PRINCIPAL COMPONENTS
pdf(paste0('figures/USE_corr_plots_PCs_', Sys.Date(), '.pdf'), height = 6, width = 8)

for (exp in c('env', 'nic', 'dis', 'mut', 'com')) {
  
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 2, byrow = F))
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), cex.axis = 1.3)
  
  for (met in c('PC1', 'PC2', 'PC3', 'PC4')) {
    
    tmp = filter(corrOutput2, experiment == exp, metric == met)
    
    plot(tmp$r, 1:nrow(tmp), pch = 18, col = tmp$color, 
         cex = 1.8, xlim = c(-1, 1), xlab = '', yaxt = 'n', ylab = '', ylim = c(0.5, 1.1*nrow(tmp)))
    abline(v = 0, col = 'black', lwd = 2)
    segments(tmp$r.L95, 1:nrow(tmp), tmp$r.U95, 1:nrow(tmp), col = tmp$color, lwd = 2.5)
    text(-1, nrow(tmp), met, cex = 1.5, adj = c(0, 1))
    
  }
  # legend panel
  plot(1, 1, type = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', bty = 'n')
  
  points(rep(0.58, nrow(tmp)), seq(0.62, 1.38, length.out = nrow(tmp)),
         pch = 18, col = tmp$color, cex = 2)
  text(rep(.61, nrow(tmp)), seq(0.62, 1.38, length.out = nrow(tmp)), 
       paste(tmp$model, '-', tmp$parameter), cex = 1.3, adj = 0)
  
  mtext(paste(experiments$phrase[experiments$experiment == exp], "experiment"), 3, outer = T, cex = 2, line = 1)
  mtext("correlation coefficient", 1, outer = T, cex = 1.5, line = 0, at = .3)
  
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), mfrow = c(4, 4))
}
dev.off()










# Plot of distribution of correlation coefficients by model
nModels = length(unique(corrOutput2$model))

if (nModels < 9) {
  
  pdf(paste0('figures/USE_r_by_model_', Sys.Date(), '.pdf'), height = 8, width = 6)
  par(mfrow = c(nModels, 1), mar = c(3, 3, 1, 1), oma = c(4, 0, 0, 0))
  
} else {
  
  pdf(paste0('figures/USE_r_by_model_', Sys.Date(), '.pdf'), height = 8, width = 10)
  par(mfrow = c(ceiling(nModels/2), 2), mar = c(3, 3, 1, 1), oma = c(4, 0, 0, 0))
  
}

for (mod in unique(corrOutput2$model)) {
  hist(corrOutput2$r[corrOutput2$model == mod], xlab = '', main = '', col = corrOutput2$color[corrOutput2$model == mod], xlim = c(-1, 1), breaks = seq(-1, 1, by = .1))
  legend("topleft", mod, bty = 'n', cex = 2)
}
mtext("Correlation coefficient", 1, cex = 2, outer = T, line = 1)
dev.off()




# Plot of distribution of correlation coefficients by experiment
pdf(paste0('figures/USE_r_by_experiment_', Sys.Date(), '.pdf'), height = 8, width = 6)
par(mfrow = c(5, 1), mar = c(3, 3, 1, 1), oma = c(4, 0, 0, 0))
for (exp in c('env', 'nic', 'dis', 'mut', 'com')) {
  hist(corrOutput2$r[corrOutput2$experiment == exp], xlab = '', main = '', col = 'gray50', 
       xlim = c(-1, 1), breaks = seq(-1, 1, by = .1))
  legend("topleft", exp, bty = 'n', cex = 2)
  abline(v = mean(corrOutput2$r[corrOutput2$experiment == exp], na.rm = T), col = 'red')
}
mtext("Correlation coefficient", 1, cex = 2, outer = T, line = 1)
dev.off()

