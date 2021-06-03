# Analyzing Uniform Sampling Experiment trees

# Load libraries
library(gsheet)
library(dplyr)
library(stringr)
library(ape)
library(geiger)
library(lessR)
library(gsheet)

source('code/pcaFunctions.r')
source('code/experiment_analysis_functions.r')
source('code/treeMetrics.r')



# Script for calculating tree metrics for trees that have not already been analyzed

# Read in existing output file
treeOutput = read.table('USE_treeOutput.txt', sep = '\t', header = T)

analyzedTrees = paste(treeOutput$model, "_", treeOutput$simID, ".tre", sep = "")

# Get list of tree files in USE experiment folder
allTreeFiles = list.files('trees/uniform_sampling_experiment')[grepl(".tre", list.files('trees/uniform_sampling_experiment'))]

treesToRun = allTreeFiles[!allTreeFiles %in% analyzedTrees]

# Calculate tree metrics for remaining trees and append
metricsForManyTrees(treefiles = treesToRun, minimumTreeSize = 5, fileOut = 'USE_treeOutput.txt', append = TRUE, 
                    treedir = 'trees/uniform_sampling_experiment')



# Function for calculating CI for spearman's rho
# --from https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation
spearman_CI <- function(x, y, alpha = 0.05){
  rs <- cor(x, y, method = "spearman", use = "complete.obs")
  n <- sum(complete.cases(x, y))
  sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
}





# For each process (env, nic, dis, mut, tim), examine the relationship between
# the strength of the parameter value associated with that process and tree metrics.

# 

# Calculate correlations between tree metrics and treatment levels for each model
corrCalcUSE = function(experiment, experimentData, modelAbbrev, cor.method = 'spearman') {
  
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
                          use = 'na.or.complete', method = cor.method)
          
          corDF$r[m + (p-1)*length(metrics)] = sign*corr$estimate
          
          if (cor.method == 'pearson') {
            
            corDF$r.L95[m + (p-1)*length(metrics)] = sign*corr$conf.int[1]
            corDF$r.U95[m + (p-1)*length(metrics)] = sign*corr$conf.int[2]
            
          } else if (cor.method == 'spearman') {
            
            corCI = spearman_CI(modelData[, metrics[m]], modelData[, experimentalParam[p]])
            corDF$r.L95[m + (p-1)*length(metrics)] = sign*corCI[1]
            corDF$r.U95[m + (p-1)*length(metrics)] = sign*corCI[2]
            
          }
          
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
                                    'mutation/speciation rate', 'competition', 'time'))

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
    
    # Check for new model-experiment combinations to run
    if (!paste(mod, e) %in% paste(corrOutput$model, corrOutput$experiment)) {

      corDF = corrCalcUSE(e, metrics, mod)
      
      if (class(corDF) == 'data.frame') {
        corrOutput = rbind(corrOutput, corDF)
      }
      
    }

  }
    
}


models = unique(paramKey$model)
modelColors = data.frame(model = models, color = colorSelection(length(models)))
modelColors$model = as.character(modelColors$model)
modelColors$color = as.character(modelColors$color)

corrOutput2 = left_join(corrOutput, modelColors, by = 'model')




# Plotting correlation coefficients by model
pdf('figures/USE_corr_plots.pdf', height = 8, width = 10)

for (exp in c('env', 'nic', 'dis', 'mut', 'com', 'tim')) {
  
  layout(matrix(c(1:15, 15), nrow = 4, byrow = T))
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), cex.axis = 1.3)
  
  for (met in unique(corrOutput2$metric)) {

    tmp = filter(corrOutput2, experiment == exp, metric == met)
    
    plot(tmp$r, 1:nrow(tmp), pch = 18, col = tmp$color, 
         cex = 3, xlim = c(-1, 1), xlab = '', yaxt = 'n', ylab = '', ylim = c(0.5, 1.1*nrow(tmp)))
    segments(tmp$r.L95, 1:nrow(tmp), tmp$r.U95, 1:nrow(tmp), col = tmp$color, lwd = 2)
    text(-1, nrow(tmp), met, cex = 1.5, adj = c(0, 1))
    abline(v = 0, col = 'black', lwd = 2)
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
         paste(tmp$model[1:firstHalf], '-', tmp$parameter[1:firstHalf]), cex = 1.3, adj = 0)

    points(rep(1.08, (numMods - firstHalf)), seq(0.62, 1.38, length.out = (numMods - firstHalf)),
           pch = 18, col = tmp$color[(firstHalf + 1):numMods], cex = 2)
    text(rep(1.12, (numMods - firstHalf)), seq(0.62, 1.38, length.out = (numMods - firstHalf)), 
         paste(tmp$model[(firstHalf + 1):numMods], '-', tmp$parameter[(firstHalf + 1):numMods]), cex = 1.3, adj = 0)
    
  }

  mtext(paste(experiments$phrase[experiments$experiment == exp], "experiment"), 3, outer = T, cex = 2, line = 1)
  mtext("correlation coefficient", 1, outer = T, cex = 1.5, line = 0, at = .3)
  
  par(mar = c(4, 4, 0, 1), oma = c(3, 0, 4, 0), mgp = c(2.2, 1, 0), mfrow = c(4, 4))
}
dev.off()


# Plot of distribution of correlation coefficients by model
pdf('figures/USE_r_by_model.pdf', height = 8, width = 6)
par(mfrow = c(length(unique(corrOutput2$model)), 1), mar = c(3, 3, 1, 1), oma = c(4, 0, 0, 0))
for (mod in unique(corrOutput2$model)) {
  hist(corrOutput2$r[corrOutput2$model == mod], xlab = '', main = '', col = corrOutput2$color[corrOutput2$model == mod], xlim = c(-1, 1), breaks = seq(-1, 1, by = .1))
  legend("topleft", mod, bty = 'n', cex = 2)
}
mtext("Correlation coefficient", 1, cex = 2, outer = T, line = 1)
dev.off()


# Plot of distribution of correlation coefficients by experiment
pdf('figures/USE_r_by_experiment.pdf', height = 8, width = 6)
par(mfrow = c(5, 1), mar = c(3, 3, 1, 1), oma = c(4, 0, 0, 0))
for (exp in c('env', 'nic', 'dis', 'mut', 'com')) {
  hist(corrOutput2$r[corrOutput2$experiment == exp], xlab = '', main = '', col = 'gray50', 
       xlim = c(-1, 1), breaks = seq(-1, 1, by = .1))
  legend("topleft", exp, bty = 'n', cex = 2)
}
mtext("Correlation coefficient", 1, cex = 2, outer = T, line = 1)
dev.off()

