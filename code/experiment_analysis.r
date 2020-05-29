# Join parameter files indicating high, medium, and low treatment levels to tree output.

library(dplyr)

source('code/pcaFunctions.r')

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
joinedOutput = joinedOutput[-1, ] %>%
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
           timStandardized == 'H' ~ 3))
         



## Plotting results

# 1. For a given experiment, plot model means for low, med, and high values

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



# Statistical analysis - 1. correlations

corrCalc = function(experiment, experimentData, mod) {
  
  if (! experiment %in% c('env', 'nic', 'dis', 'mut', 'tim')) {
    stop("'experiment' but be either 'env', 'nic', 'dis', 'mut', or 'tim'.")
  }
  
  modelData = filter(experimentData, model2 == mod)
  
  # Only calculate a correlation coefficient if there are 3 distinct levels for a given experiment
  if (sum(1:3 %in% unique(modelData[, paste(experiment, "Level", sep = "")])) == 3) {
    corDF = data.frame(model = mod, 
                        experiment = experiment,
                        r.log10S = cor(modelData$log10S, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.PD = cor(modelData$PD, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.Gamma = cor(modelData$Gamma, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.Beta = cor(modelData$Beta, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.Colless = cor(modelData$Colless, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.Sackin = cor(modelData$Sackin, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.Yule.PDA.ratio = cor(modelData$Yule.PDA.ratio, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.MRD = cor(modelData$MRD, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.VRD = cor(modelData$VRD, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.PSV = cor(modelData$PSV, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.mean.Iprime = cor(modelData$mean.Iprime, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.MPD = cor(modelData$MPD, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.VPD = cor(modelData$VPD, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'),
                        r.nLTT_stat = cor(modelData$nLTT_stat, modelData[,paste(experiment, "Level", sep = "")], use = 'complete.obs', method = 'spearman'))
      
  } else {
    corDF = data.frame(model = mod, 
                        experiment = experiment,
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
    
  }
    
  return(corDF)  
}



corrOutput = data.frame(model = NA, 
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

for (m in unique(joinedOutput$model2)) {
  
  for (e in c('env', 'nic', 'dis', 'mut', 'tim')) {
    corrs = corrCalc(experiment = e, joinedOutput, mod = m)
  
    corrOutput = rbind(corrOutput, corrs)
  }
}



# Statistical analysis 2. Mixed effects model

envMod = lmer(metric ~ envLevel + (1 + envLevel|model2), data = joinedOutput, REML = FALSE)
coef()
