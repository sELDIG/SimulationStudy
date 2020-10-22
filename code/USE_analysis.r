# Analyzing Uniform Sampling Experiment trees

# For each process (env, nic, dis, mut, tim), examine the relationship between
# the strength of the parameter value associated with that process and tree metrics.

# 

# Load libraries
library(gsheet)



# Calculate correlations between tree metrics and treatment levels for each model
corrCalc = function(experiment, experimentData, mod) {
  
  if (! experiment %in% c('env', 'nic', 'dis', 'mut', 'tim')) {
    stop("'experiment' but be either 'env', 'nic', 'dis', 'mut', or 'tim'.")
  }
  
  url <- 'https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1171496897'
  paramKey <- gsheet2tbl(url)
  
  modelParams = read.csv(paste0('trees/uniform_sampling_experiment/', model, '_USE_parameters.csv'), header = T)
  
  modelData = filter(experimentData, model2 == mod) %>%
    left_join(modelParams, by = 'simID')
  
  
  
  experimentalParam = paramKey$parameterName[paramKey$model == mod & paramKey$process == experiment]
  
    corDF = data.frame(model2 = mod, 
                       experiment = experiment,
                       r.log10S = cor(modelData$log10S, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.PD = cor(modelData$PD, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.Gamma = cor(modelData$Gamma, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.Beta = cor(modelData$Beta, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.Colless = cor(modelData$Colless, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.Sackin = cor(modelData$Sackin, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.Yule.PDA.ratio = cor(modelData$Yule.PDA.ratio, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.MRD = cor(modelData$MRD, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.VRD = cor(modelData$VRD, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.PSV = cor(modelData$PSV, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.mean.Iprime = cor(modelData$mean.Iprime, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.MPD = cor(modelData$MPD, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.VPD = cor(modelData$VPD, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'),
                       r.nLTT_stat = cor(modelData$nLTT_stat, modelData[,paste(experiment, "Level", sep = "")], use = 'na.or.complete', method = 'spearman'))
    
  } else {
    corDF = NA
  }
  
  return(corDF)  
}
