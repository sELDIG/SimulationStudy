# Functions for analyzing and summarizing simulation experiment data

# Plot tree metrics across low, medium, and high levels for a given experiment
plotExperimentResults = function(
  experiment,                       #specify 'env', 'nic', 'dis', 'mut', or 'tim'
  metrics,                          #specify dataframe with tree metrics and treatment level columns ('joinedOutput' above)
  alpha = 255                       #color transparency (255 = opaque, 1 = transparent)
) {
  
  if (!experiment %in% c('env', 'nic', 'dis', 'mut', 'tim')) {
    stop("'experiment' must be one of the following abbreviations: env, nic, dis, mut, or tim")
  }
  
  if (! experiment %in% names(metrics)) {
    stop("The 'metrics' dataframe does not have a column for that experiment.")
  }
  
  names(metrics)[names(metrics) == experiment] = 'treatment'
  
  modelColors = data.frame(model2 = unique(metrics$model2), color = colorSelection(length(unique(metrics$model2)), alpha = alpha))
  modelColors$model2 = as.character(modelColors$model2)
  modelColors$color = as.character(modelColors$color)
  
  experiments = data.frame(experiment = c('env', 'nic', 'dis', 'mut', 'tim'), 
                           phrase = c('environmental filtering', 'niche conservatism', 'disperal', 'mutation/speciation rate', 'time'))
  
  grouped = metrics %>%
    filter(treatment != "") %>%
    group_by(treatment, model2) %>%
    summarize(n = n(),
              mean_log10S = mean(log10S, na.rm = T), sd_log10S = var(log10S, na.rm = T)^.5,
              mean_PD = mean(PD, na.rm = T), sd_PD = var(PD, na.rm = T)^.5,
              mean_Gamma = mean(Gamma, na.rm = T), sd_Gamma = var(Gamma, na.rm = T)^.5,
              mean_Beta = mean(Beta, na.rm = T), sd_Beta = var(Beta, na.rm = T)^.5,
              mean_Colless = mean(Colless, na.rm = T), sd_Colless = var(Colless, na.rm = T)^.5,
              mean_Sackin = mean(Sackin, na.rm = T), sd_Sackin = var(Sackin, na.rm = T)^.5,
              mean_Yule.PDA.ratio = mean(Yule.PDA.ratio, na.rm = T), sd_Yule.PDA.ratio = var(Yule.PDA.ratio, na.rm = T)^.5,
              mean_MRD = mean(MRD, na.rm = T), sd_MRD = var(MRD, na.rm = T)^.5,
              mean_VRD = mean(VRD, na.rm = T), sd_VRD = var(VRD, na.rm = T)^.5,
              mean_PSV = mean(PSV, na.rm = T), sd_PSV = var(PSV, na.rm = T)^.5,
              mean_mean.Iprime = mean(mean.Iprime, na.rm = T), sd_mean.Iprime = var(mean.Iprime, na.rm = T)^.5,
              mean_MPD = mean(MPD, na.rm = T), sd_MPD = var(MPD, na.rm = T)^.5,
              mean_VPD = mean(VPD, na.rm = T), sd_VPD = var(VPD, na.rm = T)^.5,
              mean_nLTT_stat = mean(nLTT_stat, na.rm = T), sd_nLTT_stat = var(nLTT_stat, na.rm = T)^.5  ) %>%
    mutate(treatmentStandardized = toupper(substr(treatment, 1, 1)),
           level = case_when(
             treatmentStandardized == 'L' ~ 1,
             treatmentStandardized == 'M' ~ 2,
             treatmentStandardized == 'H' ~ 3)) %>%
    arrange(level) %>%
    left_join(modelColors, by = 'model2') %>%
    dplyr::select(model2, treatment, level, mean_log10S:sd_nLTT_stat, color)
  
  
  pdf(paste('figures/', experiment, '_results_', Sys.Date(), '.pdf', sep = ''), height = 8, width = 10)
  par(mfrow = c(4, 4), mar = c(3, 4, 1, 1), oma = c(0, 0, 3, 0), mgp = c(2.5, 1, 0), cex.lab = 1.5)
  
  for (p in 1:14) {
    plot(grouped$level, c(as.matrix(grouped[, 2*p+2])), type = 'n', xaxt = 'n', 
         xlab = '', ylab = strsplit(names(grouped)[2*p+2], "_")[[1]][2])
    mtext(c('Low', 'Med', 'High'), 1, at = 1:3, line = 1)
    
    for (m in unique(grouped$model2)) {
      points(grouped$level[grouped$model2 == m], c(as.matrix(grouped[grouped$model2 == m, 2*p+2])), 
             type = 'b', col = grouped$color[grouped$model2 == m], pch = 16, cex = 2)
      #lines()
      
    }
  }
  
  plot(1, 1, type = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', bty = 'n')
  points(rep(0.8, nrow(modelColors)), seq(0.7, 1.3, length.out = nrow(modelColors)),
         pch = 16, col = modelColors$color, cex = 3)
  text(rep(1, nrow(modelColors)), seq(0.7, 1.3, length.out = nrow(modelColors)), modelColors$model2)
  
  plot(1, 1, type = 'n', xlab = '', ylab = '', yaxt = 'n', xaxt = 'n', bty = 'n')
  
  mtext(paste("Experimental results varying", experiments$phrase[experiments$experiment == experiment]), 
        cex = 2, outer = T)
  
  dev.off()
  
  
}




# Calculate correlations between tree metrics and treatment levels for each model
corrCalc = function(experiment, experimentData, mod) {
  
  if (! experiment %in% c('env', 'nic', 'dis', 'mut', 'tim')) {
    stop("'experiment' but be either 'env', 'nic', 'dis', 'mut', or 'tim'.")
  }
  
  modelData = filter(experimentData, model2 == mod)
  
  # Only calculate a correlation coefficient if there are 3 distinct levels for a given experiment
  if (sum(1:3 %in% unique(modelData[, paste(experiment, "Level", sep = "")])) == 3) {
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
corrCalcUSE = function(experimentData, 
                       experiment, 
                       modelAbbrev, 
                       cor.method = 'spearman',
                       metrics = c('log10S', 'PD', 'Gamma', 'Beta', 'Colless', 'Sackin', 'Yule.PDA.ratio', 'MRD', 'VRD', 'PSV',
                                   'mean.Iprime', 'MPD', 'VPD', 'nLTT_stat')) {
  
  require(stringr)
  
  if (! experiment %in% c('env', 'nic', 'dis', 'mut', 'tim', 'com')) {
    stop("'experiment' but be either 'env', 'nic', 'dis', 'mut', 'com', or 'tim'.")
  }
  
  paramKey <- read.csv('experiments/uniform_sampling_experiment/simulation_parameters_key.csv', header = T)
  
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
    right_join(modelParams, by = c('simID', 'model', 'model2')) %>%
    filter(model2 == modelAbbrev)
  
  # Check for simulation output of the specified model for this experiment
  experimentalParam = paramKey$parameterName[paramKey$model == modelAbbrev & paramKey$experiment == experiment]
  
  if (length(experimentalParam) > 0) { #check that there is a parameter for this experiment
    
    
    
    corDF = data.frame(model = rep(modelAbbrev, length(metrics)*length( experimentalParam)),
                       experiment = rep(experiment, length(metrics)*length( experimentalParam)),
                       parameter = rep(experimentalParam, each = length(metrics)),
                       metric = rep(metrics, length( experimentalParam)),
                       r = rep(NA, length(metrics)*length( experimentalParam)),
                       r.L95 = rep(NA, length(metrics)*length( experimentalParam)),
                       r.U95 = rep(NA, length(metrics)*length( experimentalParam)))
    
    for (p in 1:length(experimentalParam)) {   # if two params are listed, then cycle through
      
      # Some parameters may be configured such that there is a positive correlation between the
      # strength of the process and the parameter value and vice versa. Multiply correlations through
      # by this 'sign' so that the interpretation is similar for all parameters:
      
      # A positive correlation means that the process increases in strength
      sign = paramKey$sign[paramKey$model == modelAbbrev & paramKey$parameterName == experimentalParam[p]]
      
      
      for (m in 1:length(metrics)) { # conduct a correlation with the model parameter and each tree metric
        
        if (sum(!is.na(modelData[, metrics[m]])) > 3) { # require at least 4 data points to calculate CI's
          
          corr = cor.test(modelData[, metrics[m]], modelData[, paste0(experiment, p)], 
                          use = 'na.or.complete', method = cor.method)
          
          corDF$r[m + (p-1)*length(metrics)] = sign*corr$estimate
          
          if (cor.method == 'pearson') {
            
            corDF$r.L95[m + (p-1)*length(metrics)] = sign*corr$conf.int[1]
            corDF$r.U95[m + (p-1)*length(metrics)] = sign*corr$conf.int[2]
            
          } else if (cor.method == 'spearman') {
            
            corCI = spearman_CI(modelData[, metrics[m]], modelData[, paste0(experiment, p)])
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



# Function for extracting salient features of model output from an lmer model
extractLMEoutput = function(lmeObject, expName, metricName) {
  
  lme.summ = summary(lmeObject)
  
  out.df = data.frame(experiment = rep(expName, lme.summ$ngrps + 1), 
                      metric = rep(metricName, lme.summ$ngrps + 1),
                      model = c('global', row.names(coef(lmeObject)[[1]])),
                      slope = c(coef(lme.summ)[2, 1], coef(lmeObject)[[1]][, 2]),
                      intercept = c(coef(lme.summ)[1, 1], coef(lmeObject)[[1]][, 1]),
                      ngroups = rep(lme.summ$ngrps, lme.summ$ngrps + 1),
                      totalObs = rep(lme.summ$devcomp$dims[1], lme.summ$ngrps + 1),
                      global.slope.sd = c(coef(lme.summ)[2, 2], rep(NA, lme.summ$ngrps)),
                      global.int.sd = c(coef(lme.summ)[1, 2], rep(NA, lme.summ$ngrps)))
  
  return(out.df)
  
}





# Function for aligning simulation parameters across models according to the process each parameter is associated with.
# I.e., create columns for each process (in some cases, two columns for a process because some models have two parameters
# associated with that process) in which the relevant parameter value is stored.

# Associations between parameters and processes is originally given here: 
# https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1171496897
# but has been written to 'experiments/uniform_sampling_experiment/simulation_parameters_key.csv'

# Parameters are multiplied by the sign specified in the above spreadsheet to ensure that the strength of the process
# increases with an increase in the parameter value.

alignParametersWithProcesses = function(modelAbbrev) {
  
  params = read.csv(paste("trees/uniform_sampling_experiment/", modelAbbrev, "_USE_parameters.csv", sep = ""), header = T)
  
  if ("scenario" %in% names(params)) {
    params$model2 = paste(modelAbbrev, ".", params$scenario, sep = "")
  } else {
    params$model2 = params$model
  }
  
  # Here we drop the scenario description, assuming that the process-parameter association is not scenario-dependent
  paramKey <- read.csv('experiments/uniform_sampling_experiment/simulation_parameters_key.csv', header = T) %>% 
    mutate(model = word(model, sep = "\\.")) %>%
    filter(model == modelAbbrev) %>%
    distinct()
  
  # parameter names associated with each process
  env1name = ifelse("env" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "env"][1], NA)
  env2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "env"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "env"][2], NA)
  dis1name = ifelse("dis" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "dis"][1], NA)
  dis2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "dis"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "dis"][2], NA)
  nic1name = ifelse("nic" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "nic"][1], NA)
  nic2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "nic"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "nic"][2], NA)
  mut1name = ifelse("mut" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "mut"][1], NA)
  mut2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "mut"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "mut"][2], NA)
  com1name = ifelse("com" %in% paramKey$experiment, paramKey$parameterName[paramKey$experiment == "com"][1], NA)
  com2name = ifelse(length(unique(paramKey$parameterName[paramKey$experiment == "com"])) == 2, 
                    paramKey$parameterName[paramKey$experiment == "com"][2], NA)
  
  
  outputDF = params %>%
    mutate(
      env1Name = env1name,
      env2Name = env2name,
      dis1Name = dis1name,
      dis2Name = dis2name,
      nic1Name = nic1name,
      nic2Name = nic2name,
      mut1Name = mut1name,
      mut2Name = mut2name,
      com1Name = com1name,
      com2Name = com2name
    )
    
  outputDF$env1 = ifelse(!is.na(outputDF$env1Name), outputDF[, env1name], NA)
  outputDF$env2 = ifelse(!is.na(outputDF$env2Name), outputDF[, env2name], NA)
  outputDF$dis1 = ifelse(!is.na(outputDF$dis1Name), outputDF[, dis1name], NA)
  outputDF$dis2 = ifelse(!is.na(outputDF$dis2Name), outputDF[, dis2name], NA)
  outputDF$nic1 = ifelse(!is.na(outputDF$nic1Name), outputDF[, nic1name], NA)
  outputDF$nic2 = ifelse(!is.na(outputDF$nic2Name), outputDF[, nic2name], NA)
  outputDF$mut1 = ifelse(!is.na(outputDF$mut1Name), outputDF[, mut1name], NA)
  outputDF$mut2 = ifelse(!is.na(outputDF$mut2Name), outputDF[, mut2name], NA)
  outputDF$com1 = ifelse(!is.na(outputDF$com1Name), outputDF[, com1name], NA)
  outputDF$com2 = ifelse(!is.na(outputDF$com2Name), outputDF[, com2name], NA)
  
  output = outputDF %>%
    dplyr::select(model, model2, simID, env1Name:com2)
  
  return(output)  
}
