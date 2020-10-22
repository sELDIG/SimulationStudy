# Analyzing Uniform Sampling Experiment trees

# For each process (env, nic, dis, mut, tim), examine the relationship between
# the strength of the parameter value associated with that process and tree metrics.

# 

# Load libraries
library(gsheet)
library(dplyr)
library(stringr)

source('code/pcaFunctions.r')
source('code/experiment_analysis_functions.r')

# Calculate correlations between tree metrics and treatment levels for each model
corrCalcUSE = function(experiment, experimentData, modelAbbrev) {
  
  if (! experiment %in% c('env', 'nic', 'dis', 'mut', 'tim', 'com')) {
    stop("'experiment' but be either 'env', 'nic', 'dis', 'mut', 'com', or 'tim'.")
  }
  
  url <- 'https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1171496897'
  paramKey <- gsheet2tbl(url)
  
  # split out the root model abbreviation if necessary
  mod = word(modelAbbrev, 1, sep = fixed("."))
  
  modelParams = read.csv(paste0('trees/uniform_sampling_experiment/', mod, '_USE_parameters.csv'), header = T)
  
  
  # In some cases a single simulation model can be run under different "scenarios" that should
  # effectively be analyzed separately. In this case, a "scenario" column should be added to
  # the [model]_USE_parameters.csv file with a code for distinguishing them.

  if ('scenario' %in% names(modelParams)) {
    modelParams$model2 = paste(modelParams$model, modelParams$scenario, sep = '.')  
  } else {
    modelParams$model2 = modelParams$model
  }
  
  modelData = experimentData %>%
    right_join(modelParams, by = c('simID', 'model')) %>%
    filter(model2 == modelAbbrev)
  
  experimentalParam = paramKey$parameterName[paramKey$model == modelAbbrev & paramKey$experiment == experiment]
  
  if (length(experimentalParam) > 0) { #check that there is a parameter for this experiment

    metrics = c('log10S', 'PD', 'Gamma', 'Beta', 'Colless', 'Sackin', 'Yule.PDA.ratio', 'MRD', 'VRD', 'PSV',
                'mean.Iprime', 'MPD', 'VPD', 'nLTT_stat')
    
    corDF = data.frame(model = rep(modelAbbrev, length(metrics)*length( experimentalParam)),
                       experiment = rep(experiment, length(metrics)*length( experimentalParam)),
                       parameter = rep(experimentalParam, each = length(metrics)),
                       metric = rep(metrics, length( experimentalParam)),
                       r = rep(NA, length(metrics)*length( experimentalParam)),
                       r.L95 = rep(NA, length(metrics)*length( experimentalParam)),
                       r.U95 = rep(NA, length(metrics)*length( experimentalParam)))
    
    for (p in 1:length(experimentalParam)) {   # if multiple params are listed, then cycle through
      
      # Some parameters may be configured such that there is a positive correlation between the
      # strength of the process and the parameter value and vice versa. Multiply correlations through
      # by this 'sign' so that the interpretation is similar for all parameters:
      
      # A positive correlation means that the process increases in strength
      sign = paramKey$sign[paramKey$model == modelAbbrev & paramKey$parameterName == experimentalParam[p]]
      

      for (m in 1:length(metrics)) { # conduct a correlation with the model parameter and each tree metric
        
        if (sum(!is.na(modelData[, metrics[m]])) > 3) { # require at least 4 data points to calculate CI's
          
          corr = cor.test(modelData[, metrics[m]], modelData[, experimentalParam[p]], 
                          use = 'na.or.complete', method = 'pearson')
          
          corDF$r[m + (p-1)*length(metrics)] = sign*corr$estimate
          corDF$r.L95[m + (p-1)*length(metrics)] = sign*corr$conf.int[1]
          corDF$r.U95[m + (p-1)*length(metrics)] = sign*corr$conf.int[2]
          
        } else {
          
          corDF$r[m + (p-1)*length(metrics)] = NA
          corDF$r.L95[m + (p-1)*length(metrics)] = NA
          corDF$r.U95[m + (p-1)*length(metrics)] = NA
          
        }
        
      } # end metric loop
      
    } # end parameter loop
        
  } else { # if there is no parameter associated with this experiment, then NA
    
    corDF = NA
    
  }

  return(corDF)  
}



#########################################################################################
# Generate correlations between all parameters relevant to experiments and tree metrics

experiments = data.frame(experiment = c('env', 'nic', 'dis', 'mut', 'com', 'tim'), 
                         phrase = c('environmental filtering', 'niche conservatism', 'dispersal', 
                                    'competition', 'mutation/speciation rate', 'time'))

url <- 'https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1171496897'
paramKey <- gsheet2tbl(url)

# Read in output and join to parameter files
metrics = read.table('USE_treeOutput.txt', header = T, sep = '\t')

corrOutput = data.frame(model = character(),
                        experiment = character(),
                        parameter = character(),
                        metric = character(),
                        r = double(),
                        r.L95 = double(),
                        r.U95 = double())

for (e in experiments$experiment) {

  relevantModels = filter(paramKey, experiment == e) %>%
    distinct(model) %>% unlist()
  
  for (mod in relevantModels) {
    
    corDF = corrCalcUSE(e, metrics, mod)
    
    if (class(corDF) == 'data.frame') {
      corrOutput = rbind(corrOutput, corDF)
    }
    
  }
    
}


models = unique(paramKey$model)
modelColors = data.frame(model = models, color = colorSelection(length(models)))
modelColors$model = as.character(modelColors$model)
modelColors$color = as.character(modelColors$color)

corrOutput = left_join(corrOutput, modelColors, by = 'model')




# Plotting correlation coefficients by model
pdf('figures/USE_corr_plots.pdf', height = 8, width = 10)
par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), mfrow = c(4, 4), cex.axis = 1.3)

for (exp in c('env', 'nic', 'dis', 'mut', 'tim', 'com')) {
  
  for (met in unique(corrOutput$metric)) {
    tmp = filter(corrOutput, experiment == exp, metric == met)
    plot(tmp$r, 1:nrow(tmp), pch = 18, col = tmp$color, 
         cex = 2, xlim = c(-1, 1), xlab = '', yaxt = 'n', ylab = '')
    segments(tmp$r.L95, 1:nrow(tmp), tmp$r.U95, 1:nrow(tmp), col = tmp$color, lwd = 2)
    text(-1, nrow(tmp), met, cex = 1.5, adj = c(0, 1))
    abline(v = 0, col = 'black', lwd = 2)
  }
  # legend panel
  plot(1, 1, type = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', bty = 'n')
  points(rep(0.8, nrow(tmp)), seq(0.6, 1.4, length.out = nrow(tmp)),
         pch = 18, col = tmp$color, cex = 2)
  text(rep(1, nrow(tmp)), seq(0.6, 1.4, length.out = nrow(tmp)), tmp$model, cex = 1.5)
  
  
  mtext(paste(experiments$phrase[experiments$experiment == exp], "experiment"), 3, outer = T, cex = 2, line = 1)
  mtext("correlation coefficient", 1, outer = T, cex = 1.5, line = 0, at = .3)
  
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), mfrow = c(4, 4))
}
dev.off()




