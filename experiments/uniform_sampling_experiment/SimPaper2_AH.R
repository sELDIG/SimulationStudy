# This script 
# 1) examines simulated trees and metrics and relates them to simulation parameter values.
# 2) identifies which process-parameters are most correlated with which tree metrics
# 3) conducts random forest models for describing which metrics are most predictive of parameter values
# 4) inverts the random forest models to infer estimated parameter values given tree shape metrics for several pairs of sister clades

library(RColorBrewer)
library(dplyr)

# Read in simulated tree data
simuData <- read.csv("./experiments/uniform_sampling_experiment/process_parameter_values_and_tree_metrics_sign_corrected.csv", stringsAsFactors = T)
simuDataKey <- read.csv("./experiments/uniform_sampling_experiment/simulation_parameters_key.csv")


head(simuDataKey)

# All tree metrics
statistics = colnames(simuData)[25:43]
statisticsIndices = 25:43

# Only the tree metrics not strongly correlated with richness
#   --excluding log10S, PD, Colless, Sackin, Yule.PDA.ratio
statistics2 = colnames(simuData)[c(26, 28, 29, 33:43)]
statisticsIndices2 = c(26, 28, 29, 33:43)

# Indices for the simulation processes (env, dis, mut, nic, com)
predictorsIndices = c(14,16,18,20,22)
predictors = colnames(simuData)[c(14,16,18,20,22)]
numModelsPerPredictor = c(5, 7, 4, 7, 6)

signMat = matrix(nrow = length(statistics), ncol = length(predictors))
rownames(signMat) = statistics
colnames(signMat) = predictors
R2Mat = disAg = disAg2 = signMat

models = levels(simuData$model)


for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    R2 = sig = rep(NA, length(models))
    for(k in 1:length(models)){
      tmp = simuData[simuData$model == models[k],]
      x = scale(tmp[,statisticsIndices[i]])
      y = tmp[,predictorsIndices[j]]
      if(! all(is.na(y))){
      
      corRes = cor.test(x,y, method = "spearman")
        
      # fit <- summary(lm(y ~ x))
      # R2[k] = fit$r.squared
      # sig[k] = ifelse(fit$coefficients[2,4] <0.05, sign(fit$coefficients[2,1]), 0)
      R2[k] = corRes$estimate
      sig[k] = sign(corRes$estimate)
      } else{
        R2[k] = NA
        sig[k] = NA
      }
    }
    R2Mat[i,j] = mean(R2, na.rm = T)
    signMat[i,j] = sum(sig, na.rm = T) / sum(!is.na(sig))
    disAg[i,j] = abs(sum(sig == 1, na.rm = T) - sum(sig == -1, na.rm = T)) / abs(sum(sig %in% c(1,-1), na.rm = T))
    disAg2[i,j] = max(sum(sig == 1, na.rm = T), sum(sig == -1, na.rm = T)) / abs(sum(sig %in% c(1,-1), na.rm = T))
  }
}



corOut = data.frame(predictor = NA, statistic = NA, model = NA, r = NA)

index = 1 
for(i in 1:length(predictors)){
  for(j in 1:length(statistics)){
    for(k in 1:length(models)){
      tmp = simuData[simuData$model == models[k],]
      x = scale(tmp[,statisticsIndices[j]])
      y = tmp[,predictorsIndices[i]]
      if(! all(is.na(y))){
        
        corRes = cor.test(x,y, method = "spearman")
        
        corTmp = data.frame(predictor = predictors[i], statistic = statistics[j], model = models[k], r = round(corRes$estimate,2))
      } else{
        corTmp = data.frame(predictor = predictors[i], statistic = statistics[j], model = models[k], r = NA)
      }
      corOut = rbind(corOut, corTmp)
      index = index + 1
    }
  }
}

corrs = corOut %>%
  pivot_wider(id_cols = predictor:statistic, names_from = model, values_from = r) %>%
  slice(-1) %>%
  dplyr::select(predictor, statistic, ca:xe) %>%
  mutate(aveR = as.vector(R2Mat),
         aveSign = as.vector(signMat),
         disAg1 = as.vector(disAg),
         disAg2 = as.vector(disAg2)) %>%
  arrange(predictor, desc(disAg1))

 

image.real <- function(mat, xCol = c("blue", "white", "white", "red"), 
                       range = c(-1,1), x.labels = rownames(mat), y.labels = colnames(mat)) { 
  mat <- t(mat)[,nrow(mat):1]
  fields::image.plot(mat, axes = FALSE, zlim = range, 
                     col = colorRampPalette(xCol)(30))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = x.labels)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = y.labels, las = 2)
  box() 
}


# 2-panel figure
par(mfrow = c(1, 2), mar = c(3,11,3,3))

# Values: % agreement, Color: average sign
image.real(signMat, x.labels = c('env', 'dis', 'nic', 'mut', 'com')) 
for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices)-1), 
         labels = 100*round(disAg2[i,j], digits = 2))
  }
}
mtext(paste0(rep("(", 5), numModelsPerPredictor, rep(")", 5)), 1, at = seq(0, 1, length = ncol(signMat)), line = 2)
title(main = "% agreement in effect direction")

# Values: average correlation coefficient, Color: average sign
image.real(signMat, x.labels = c('env', 'dis', 'nic', 'mut', 'com')) 
for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices)-1), 
         labels = round(R2Mat[i,j], digits = 2))
  }
}
mtext(paste0(rep("(", 5), numModelsPerPredictor, rep(")", 5)), 1, at = seq(0, 1, length = ncol(signMat)), line = 2)
title(main = "Average correlation coefficient")


image.real(R2Mat, range = c(-1,1), xCol = c("beige", "firebrick2")) 
title(main = "R2 between parameter and tree metric")

for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices)-1), 
         labels = round(R2Mat[i,j], digits = 2))
  }
}


image.real(disAg2) 
for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices)-1), 
         labels = round(disAg[i,j], digits = 2))
  }
}


########################################################################################
# Random Forests for predicting parameter values from tree metrics
library(ranger)

# Preliminary analysis using all tree metrics

predictability = array(dim = c(length(statistics), length(predictors), length(models)))
predictabilityR2 = array(dim = c(length(predictors), length(models)))

for(j in 1:length(predictorsIndices)){
  for(k in 1:length(models)){
    tmp = simuData[simuData$model == models[k],c(predictorsIndices[j], statisticsIndices)]
    tmp = tmp[complete.cases(tmp),]
    if(nrow(tmp) > 0){
      form = as.formula(paste(predictors[j], "~ ."))
      x = ranger(form, data = tmp, importance = 'permutation',
                 scale.permutation.importance = T)
      predictability[,j,k] = importance(x)     
      predictabilityR2[j,k] = x$r.squared
    } 
  }
}


rownames(predictabilityR2) = predictors
colnames(predictabilityR2) = models

image.real(predictabilityR2, range = c(0,1), xCol = c("beige", "firebrick2"))
title(main = "Random forest predictability [R2]")


predMean = apply(predictability, c(1,2), mean, na.rm = T)
predSD = apply(predictability, c(1,2), sd, na.rm = T)

rownames(predMean) = statistics
colnames(predMean) = predictors

image.real(predMean, range = NULL, xCol = c("beige", "firebrick2"))
title(main = "mean variable importance for predicting across models")

for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices)-1), 
         labels = round(predMean[i,j] / predSD[i,j], digits = 2))
  }
}



# Analysis using only the tree metrics that are not strongly correlated with richness

predictability2 = array(dim = c(length(statistics2), length(predictors), length(models)))
predictability2R2 = array(dim = c(length(predictors), length(models)))

for(j in 1:length(predictorsIndices)){
  for(k in 1:length(models)){
    tmp = simuData[simuData$model == models[k],c(predictorsIndices[j], statisticsIndices2)]
    tmp = tmp[complete.cases(tmp),]
    if(nrow(tmp) > 0){
      form = as.formula(paste(predictors[j], "~ ."))
      x = ranger(form, data = tmp, importance = 'permutation',
                 scale.permutation.importance = T)
      predictability2[,j,k] = importance(x)     
      predictability2R2[j,k] = x$r.squared
    } 
  }
}


rownames(predictability2R2) = predictors
colnames(predictability2R2) = models

image.real(predictability2R2, range = c(0,1), xCol = c("beige", "firebrick2"))
title(main = "Random forest [R2] (non-richness metrics)")


predMean2 = apply(predictability2, c(1,2), mean, na.rm = T)
predSD2 = apply(predictability2, c(1,2), sd, na.rm = T)

rownames(predMean2) = statistics2
colnames(predMean2) = predictors

image.real(predMean2, range = NULL, xCol = c("beige", "firebrick2"))
title(main = "mean variable importance (non-richness metrics)")

for(i in 1:length(statisticsIndices2)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices2)-1), 
         labels = round(predMean2[i,j] / predSD2[i,j], digits = 2))
  }
}













#######################################################################################################
# Empirical sister clades
empiricalSisterCladeFiles = list.files('trees/empirical')[grepl('pair', list.files('trees/empirical'))]

# run a version of metricsForManyTrees() to get tree metric output (but that function assumes diff model ID naming convention,
# so had to manually tweak)
empiricalMetrics = read.table('empiricalSisterClades_treeOutput.txt', header = T, sep = '\t')
empiricalMetrics$clade = rep(c(1,2), 6)
empMetrics = empiricalMetrics[, 4:22]

predictability = array(dim = c(length(statistics), length(predictors), length(models)))
predictabilityR2 = array(dim = c(length(predictors), length(models)))


# Random forests for each model, saving predictions of parameter values from random forests for each empirical tree (based on all tree metrics)
predictedValues = data.frame(tree = character(),
                             model = character(), 
                             predictor = character(),
                             value = numeric())

for(j in 1:length(predictorsIndices)){
  for(k in 1:length(models)){
    tmp = simuData[simuData$model == models[k],c(predictorsIndices[j], statisticsIndices)]
    tmp = tmp[complete.cases(tmp),]
    if(nrow(tmp) > 0){
      form = as.formula(paste(predictors[j], "~ ."))
      x = ranger(form, data = tmp, importance = 'permutation',
                 scale.permutation.importance = T)
      predictions = predict(x, data = empMetrics)
      
      tmpPredictions = data.frame(tree = word(empiricalMetrics$tree, 1, sep = fixed('.')),
                                  model = rep(models[k], nrow(empMetrics)),
                                  predictor = rep(predictors[j], nrow(empMetrics)),
                                  value = predictions$predictions)
      
      predictedValues = rbind(predictedValues, tmpPredictions)

    } 
  }
}


# Need to standardize the predicted parameter values relative to the range of values examined within a given model.  
#   scaled = (x - min) / (max - min)          so first need to identify min and max values for each model-parameter combo:

# (OR DO QUANTILES MAKE MORE SENSE?)

modelParamRange = simuData %>%
  group_by(model) %>%
  summarize(env1Min = min(env1, na.rm = T),
            env1Max = max(env1, na.rm = T),
            com1Min = min(com1, na.rm = T),
            com1Max = max(com1, na.rm = T),
            dis1Min = min(dis1, na.rm = T),
            dis1Max = max(dis1, na.rm = T),
            nic1Min = min(nic1, na.rm = T),
            nic1Max = max(nic1, na.rm = T),
            mut1Min = min(mut1, na.rm = T),
            mut1Max = max(mut1, na.rm = T))

# Note that for some model/parameter combinations, the modeler varied parameter values uniformly on a log scale,
# while for others values varied uniformly on an arithmetic scale.

# I examined the skewness values of all sets of values, and determined that each of the following combinations
# had a skewness of 0.7 or greater and merited being log-transformed prior to scaling. 
# (Otherwise, the scaled values, e.g. in model 'hs' for dis1 or mut1 are close to 0)

# ca: mut1, com1
# hs: dis1, mut1
# xe: dis1, mut1

modelParamRange$com1Min[modelParamRange$model == 'ca'] = log10(modelParamRange$com1Min[modelParamRange$model == 'ca'])
modelParamRange$com1Max[modelParamRange$model == 'ca'] = log10(modelParamRange$com1Max[modelParamRange$model == 'ca'])
modelParamRange$dis1Min[modelParamRange$model %in% c('hs', 'xe')] = log10(modelParamRange$dis1Min[modelParamRange$model %in% c('hs', 'xe')])
modelParamRange$dis1Max[modelParamRange$model %in% c('hs', 'xe')] = log10(modelParamRange$dis1Max[modelParamRange$model %in% c('hs', 'xe')])
modelParamRange$mut1Min[modelParamRange$model %in% c('hs', 'xe', 'ca')] = log10(modelParamRange$mut1Min[modelParamRange$model %in% c('hs', 'xe', 'ca')])
modelParamRange$mut1Max[modelParamRange$model %in% c('hs', 'xe', 'ca')] = log10(modelParamRange$mut1Max[modelParamRange$model %in% c('hs', 'xe', 'ca')])

# Also need to log-transform these in predictedValues

predictedValues$value[predictedValues$model == 'ca' & predictedValues$predictor %in% c('mut1', 'com1')] = 
  log10(predictedValues$value[predictedValues$model == 'ca' & predictedValues$predictor %in% c('mut1', 'com1')])

predictedValues$value[predictedValues$model %in% c('hs', 'xe') & predictedValues$predictor %in% c('mut1', 'dis1')] = 
  log10(predictedValues$value[predictedValues$model %in% c('hs', 'xe') & predictedValues$predictor %in% c('mut1', 'dis1')])



# Rearrange
modelParamMinMax = data.frame(model = rep(modelParamRange$model, 5), 
                              predictor = rep(c('env1', 'com1', 'dis1', 'nic1', 'mut1'), each = 8),
                              minVal = c(modelParamRange$env1Min, modelParamRange$com1Min, modelParamRange$dis1Min, 
                                         modelParamRange$nic1Min, modelParamRange$mut1Min),
                              maxVal = c(modelParamRange$env1Max, modelParamRange$com1Max, modelParamRange$dis1Max, 
                                         modelParamRange$nic1Max, modelParamRange$mut1Max))



# Function for rescaling a value between the min and max such that the scaled value varies between 0 and 1.
# If either the min or max are NA, Inf, or -Inf, returns NA
scaleFunction = function(value, min, max) {
  
  if (!min %in% c(-Inf, NA, Inf) & !max %in% c(-Inf, NA, Inf)) {
    
    scaledValue = (value - min) / (max - min)
  
  } else {
  
    scaledValue = NA
  
  }
  return(scaledValue)
}



modelColors = data.frame(model = unique(predictedValues$model),
                         color = as.character(colorSelection(length(unique(predictedValues$model)))))


scaledPredictions = predictedValues %>%
  left_join(modelParamMinMax, by = c('model', 'predictor')) %>%
  mutate(pair = str_extract(tree, "[1-9]"),
         clade = rep(c(1, 2), 174))

scaledPredictions$scaledValue = apply(scaledPredictions[, c('value', 'minVal', 'maxVal')], 1, 
                             function(x) scaleFunction(value = x[1], min = x[2], max = x[3]))

scaledPredictions = left_join(scaledPredictions, modelColors, by = 'model')




# Then can plot distribution of scaled scores between sister clades (connecting points within models) for all processes
pdf('figures/empirical_sisterclade_comparisons.pdf', height = 8, width = 10)

pairs = unique(scaledPredictions$pair)
processes = unique(scaledPredictions$predictor)

for (pro in processes) {
  
  for (p in pairs) {
    
    tmp = filter(scaledPredictions, pair == p, predictor == pro)
    
    if (p == 1) {
      
      par(mfrow = c(2, 3), mar = c(3, 3, 3, 0), oma = c(0, 3, 3, 0))

    } 
    plot(c(0.5,2.5), range(tmp$scaledValue), type = 'n', xaxt = 'n', xlab = "", ylab = "", las = 1,
         main = word(tmp$tree[2], 2, sep = "_"))
    axis(1, at = 1:2, labels = c('clade 1', 'clade 2'))
    
    for (m in unique(tmp$model)) {
      
      points(tmp$clade[tmp$model==m], tmp$scaledValue[tmp$model==m], type = 'b', col = tmp$color[tmp$model==m], lwd = 3)
      text(2.1, tmp$scaledValue[tmp$model==m & tmp$clade == 2], m, adj = 0)
    }
  }

  mtext(pro, side = 3, outer = TRUE, cex = 2)
  mtext("Scaled parameter", 2, outer = TRUE, cex = 1.5)
}
dev.off()


########################################################################################
# A few example figures of a select clade and process, along with the tree metrics 
# most associated with that process based on earlier correlation analysis


# Need to pull out min and max values of each tree metric across simulations for scaling purposes
processDFmetrics = read.csv('experiments/uniform_sampling_experiment/process_parameter_values_and_tree_metrics.csv')

metricsMinMax = data.frame(metric = names(processDFmetrics)[25:43],
                           minVal = apply(processDFmetrics[, 25:43], 2, function(x) min(x, na.rm = T)),
                           maxVal = apply(processDFmetrics[, 25:43], 2, function(x) max(x, na.rm = T)))


modelParamRange = simuData %>%
  group_by(model) %>%
  summarize(env1Min = min(env1, na.rm = T),
            env1Max = max(env1, na.rm = T),
            com1Min = min(com1, na.rm = T),
            com1Max = max(com1, na.rm = T),
            dis1Min = min(dis1, na.rm = T),
            dis1Max = max(dis1, na.rm = T),
            nic1Min = min(nic1, na.rm = T),
            nic1Max = max(nic1, na.rm = T),
            mut1Min = min(mut1, na.rm = T),
            mut1Max = max(mut1, na.rm = T))






# Cetartiodactyla and dispersal
clade1 = read.tree('trees/empirical/pair2_cetartiodactyla sister clade carnivora perrisodactyla.txt')
clade2 = read.tree('trees/empirical/pair2_cetartiodactyla.txt')

pdf('figures/cetartiodactyla_dispersal.pdf', height = 8, width = 10)
par(mfcol = c(2, 2), mar = c(3, 4, 1, 1), oma = c(1, 0, 3, 0))
plot(clade1, show.tip.label = F, main = "Carnivora/Perrisodactyla", cex.main = 1.5)
plot(clade2, show.tip.label = F, main = "Cetartiodactyla", cex.main = 1.5)

# Metrics associated with dispersal are VRD, mean.Iprime, (Colless, Sackin, Yule, MRD, MGL_principal_eigenvalue)

dispMetrics = c('VRD', 'mean.Iprime', 'Beta', 'Colless', 'Sackin', 'Yule.PDA.ratio', 'MRD', 'MGL_principal_eigenvalue')


# Get scaled values of tree metrics relative to the distribution of metrics from all simulated trees
# In this case it seemed easier/more straightforward to calculate quantiles
# (in part because the distribution for many of these simulated metrics is highly skewed)
scaledMetrics = data.frame(metric = dispMetrics, 
                           clade1 = rep(NA, length(dispMetrics)), 
                           clade2 = rep(NA, length(dispMetrics)),
                           metricColor = brewer.pal(length(dispMetrics), 'Blues'))

for (d in dispMetrics) {
  
  scaledMetrics$clade1[scaledMetrics$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 2 & empiricalMetrics$clade == 1, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)
    
  scaledMetrics$clade2[scaledMetrics$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 2 & empiricalMetrics$clade == 2, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)

}



plot(c(0.5,2.5), range(c(scaledMetrics$clade1, scaledMetrics$clade2), na.rm = T), type = 'n', xaxt = 'n', xlab = "", ylab = "Scaled tree metric", las = 1, cex.lab = 1.3)
axis(1, at = 1:2, labels = c('Carnivora/\nPerrisodactyla', 'Cetartiodactyla'), tck = -.01)

for (m in scaledMetrics$metric) {
  
  points(x = c(1, 2), y = scaledMetrics[scaledMetrics$metric == m, 2:3], 
         type = 'b', col = scaledMetrics$metricColor[scaledMetrics$metric == m], lwd = 3)
  
  text(2.1, scaledMetrics$clade2[scaledMetrics$metric == m], m, adj = 0)
  
}
text(.65, max(scaledMetrics$clade1, na.rm = T), "Metrics", cex = 1.5)

# Inferred dispersal parameters for different models
dispersalParams = filter(scaledPredictions, pair == 2, predictor == 'dis1')

plot(c(0.5,2.5), range(dispersalParams$scaledValue), type = 'n', xaxt = 'n', xlab = "", ylab = "low <--    Inferred dispersal    --> high", las = 1, cex.lab = 1.3)
axis(1, at = 1:2, labels = c('Carnivora/\nPerrisodactyla', 'Cetartiodactyla'), tck = -.01)

for (m in unique(dispersalParams$model)) {
  
  points(dispersalParams$clade[dispersalParams$model==m], dispersalParams$scaledValue[dispersalParams$model==m], 
         type = 'b', col = dispersalParams$color[dispersalParams$model==m], lwd = 3)
  
  text(2.1, dispersalParams$scaledValue[dispersalParams$model==m & dispersalParams$clade == 2], m, adj = 0)
  
}
text(.65, max(dispersalParams$scaledValue, na.rm = T), "Models", cex = 1.5)

dev.off()



###############################################################################################
# Plots based on random forests using tree metrics NOT STRONGLY CORRELATED WITH RICHNESS

# Random forests for each model, saving predictions of parameter values from random forests for each empirical tree (based on tree metrics NOT STRONGLY CORRELATED WITH RICHNESS)
predictedValues2 = data.frame(tree = character(),
                              model = character(), 
                              predictor = character(),
                              value = numeric())
empMetrics2 = empiricalMetrics[, c(5, 7, 8, 12:23)]

predictability2 = array(dim = c(length(statistics2), length(predictors), length(models)))
predictability2R2 = array(dim = c(length(predictors), length(models)))


for(j in 1:length(predictorsIndices)){
  for(k in 1:length(models)){
    tmp = simuData[simuData$model == models[k],c(predictorsIndices[j], statisticsIndices2)]
    tmp = tmp[complete.cases(tmp),]
    if(nrow(tmp) > 0){
      form = as.formula(paste(predictors[j], "~ ."))
      x = ranger(form, data = tmp, importance = 'permutation',
                 scale.permutation.importance = T)
      predictions = predict(x, data = empMetrics2)
      
      tmpPredictions = data.frame(tree = word(empiricalMetrics$tree, 1, sep = fixed('.')),
                                  model = rep(models[k], nrow(empMetrics2)),
                                  predictor = rep(predictors[j], nrow(empMetrics2)),
                                  value = predictions$predictions)
      
      predictedValues2 = rbind(predictedValues2, tmpPredictions)
      
    } 
  }
}


# Need to log-transform these in predictedValues
predictedValues2$value[predictedValues2$model == 'ca' & predictedValues2$predictor %in% c('mut1', 'com1')] = 
  log10(predictedValues2$value[predictedValues2$model == 'ca' & predictedValues2$predictor %in% c('mut1', 'com1')])

predictedValues2$value[predictedValues2$model %in% c('hs', 'xe') & predictedValues2$predictor %in% c('mut1', 'dis1')] = 
  log10(predictedValues2$value[predictedValues2$model %in% c('hs', 'xe') & predictedValues2$predictor %in% c('mut1', 'dis1')])


scaledPredictions2 = predictedValues2 %>%
  left_join(modelParamMinMax, by = c('model', 'predictor')) %>%
  mutate(pair = str_extract(tree, "[1-9]"),
         clade = rep(c(1, 2), 174))

scaledPredictions2$scaledValue = apply(scaledPredictions2[, c('value', 'minVal', 'maxVal')], 1, 
                             function(x) scaleFunction(value = x[1], min = x[2], max = x[3]))

scaledPredictions2 = left_join(scaledPredictions2, modelColors, by = 'model')

 

# Then can plot distribution of scaled scores between sister clades (connecting points within models) for all processes
pdf('figures/empirical_sisterclade_comparisons_non-S-correlated_metrics.pdf', height = 8, width = 10)

pairs = unique(scaledPredictions2$pair)
processes = unique(scaledPredictions2$predictor)

for (pro in processes) {
  
  for (p in pairs) {
    
    tmp = filter(scaledPredictions2, pair == p, predictor == pro)
    
    if (p == 1) {
      
      par(mfrow = c(2, 3), mar = c(3, 3, 3, 0), oma = c(0, 3, 3, 0))
      
    } 
    plot(c(0.5,2.5), range(tmp$scaledValue), type = 'n', xaxt = 'n', xlab = "", ylab = "", las = 1,
         main = word(tmp$tree[2], 2, sep = "_"))
    axis(1, at = 1:2, labels = c('clade 1', 'clade 2'))
    
    for (m in unique(tmp$model)) {
      
      points(tmp$clade[tmp$model==m], tmp$scaledValue[tmp$model==m], type = 'b', col = tmp$color[tmp$model==m], lwd = 3)
      text(2.1, tmp$scaledValue[tmp$model==m & tmp$clade == 2], m, adj = 0)
    }
  }
  
  mtext(pro, side = 3, outer = TRUE, cex = 2)
  mtext("Scaled parameter", 2, outer = TRUE, cex = 1.5)
}
dev.off()


########################################################################################
# A few example figures of a select clade and process, along with the tree metrics 
# most associated with that process based on earlier correlation analysis


# Cetartiodactyla and dispersal
clade1 = read.tree('trees/empirical/pair2_cetartiodactyla sister clade carnivora perrisodactyla.txt')
clade2 = read.tree('trees/empirical/pair2_cetartiodactyla.txt')

pdf('figures/cetartiodactyla_dispersal_non-S-correlated_metrics.pdf', height = 8, width = 10)
par(mfcol = c(2, 2), mar = c(3, 4, 1, 1), oma = c(1, 0, 3, 0))
plot(clade1, show.tip.label = F, main = "Carnivora/Perrisodactyla", cex.main = 1.5)
plot(clade2, show.tip.label = F, main = "Cetartiodactyla", cex.main = 1.5)

# Metrics associated with dispersal are VRD, mean.Iprime, MRD, MGL_principal_eigenvalue (excluding those correlated with S)

dispMetrics2 = c('VRD', 'Beta', 'mean.Iprime', 'MRD', 'MGL_principal_eigenvalue')


# Get scaled values of tree metrics relative to the distribution of metrics from all simulated trees
# In this case it seemed easier/more straightforward to calculate quantiles
# (in part because the distribution for many of these simulated metrics is highly skewed)
scaledMetrics2 = data.frame(metric = dispMetrics2, 
                           clade1 = rep(NA, length(dispMetrics2)), 
                           clade2 = rep(NA, length(dispMetrics2)),
                           metricColor = brewer.pal(length(dispMetrics2), 'Blues'))

for (d in dispMetrics2) {
  
  scaledMetrics2$clade1[scaledMetrics2$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 2 & empiricalMetrics$clade == 1, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)
  
  scaledMetrics2$clade2[scaledMetrics2$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 2 & empiricalMetrics$clade == 2, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)
  
}



plot(c(0.5,2.5), range(c(scaledMetrics2$clade1, scaledMetrics2$clade2), na.rm = T), type = 'n', xaxt = 'n', xlab = "", ylab = "Scaled tree metric", las = 1, cex.lab = 1.3)
axis(1, at = 1:2, labels = c('Carnivora/\nPerrisodactyla', 'Cetartiodactyla'), tck = -.01)

for (m in scaledMetrics2$metric) {
  
  points(x = c(1, 2), y = scaledMetrics2[scaledMetrics2$metric == m, 2:3], 
         type = 'b', col = scaledMetrics2$metricColor[scaledMetrics2$metric == m], lwd = 3)
  
  text(2.1, scaledMetrics2$clade2[scaledMetrics2$metric == m], m, adj = 0)
  
}
text(.65, max(scaledMetrics2$clade1, na.rm = T), "Metrics", cex = 1.5)

# Inferred dispersal parameters for different models
dispersalParams2 = filter(scaledPredictions2, pair == 2, predictor == 'dis1')

plot(c(0.5,2.5), range(dispersalParams2$scaledValue), type = 'n', xaxt = 'n', xlab = "", ylab = "low <--    Inferred dispersal    --> high", las = 1, cex.lab = 1.3)
axis(1, at = 1:2, labels = c('Carnivora/\nPerrisodactyla', 'Cetartiodactyla'), tck = -.01)

for (m in unique(dispersalParams2$model)) {
  
  points(dispersalParams2$clade[dispersalParams2$model==m], dispersalParams2$scaledValue[dispersalParams2$model==m], 
         type = 'b', col = dispersalParams2$color[dispersalParams2$model==m], lwd = 3)
  
  text(2.1, dispersalParams2$scaledValue[dispersalParams2$model==m & dispersalParams2$clade == 2], m, adj = 0)
  
}
text(.65, max(dispersalParams2$scaledValue, na.rm = T), "Models", cex = 1.5)

dev.off()



########################################################################################
# Example 2: Hummingbirds vs sisters and environmental filtering

clade1 = read.tree('trees/empirical/pair4_humminbirds sister clade swifts.txt')
clade2 = read.tree('trees/empirical/pair4_hummingbirds.txt')

pdf('figures/hummingbirds_env_filtering_non-S-correlated_metrics.pdf', height = 8, width = 10)
par(mfcol = c(2, 2), mar = c(3, 4, 3, 1), oma = c(1, 0, 3, 0))
plot(clade1, show.tip.label = F, main = "", cex.main = 1.5)
mtext("A) Swifts", 3, adj = 0, cex = 1.3, bty='n', line = 1)
plot(clade2, show.tip.label = F, main = "", cex.main = 1.5)
mtext("B) Hummingbirds", 3, adj = 0, cex = 1.3, bty='n', line = 1)

# Metrics associated with environmental filtering, 80% agreement (excluding those correlated with S)

envMetrics = c('MGL_peakedness', 'Gamma', 'Beta')


# Get scaled values of tree metrics relative to the distribution of metrics from all simulated trees
# In this case it seemed easier/more straightforward to calculate quantiles
# (in part because the distribution for many of these simulated metrics is highly skewed)
scaledMetrics3 = data.frame(metric = envMetrics, 
                            clade1 = rep(NA, length(envMetrics)), 
                            clade2 = rep(NA, length(envMetrics)),
                            metricColor = brewer.pal(length(envMetrics), 'Blues'))

for (d in envMetrics) {
  
  scaledMetrics3$clade1[scaledMetrics3$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 4 & empiricalMetrics$clade == 1, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)
  
  scaledMetrics3$clade2[scaledMetrics3$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 4 & empiricalMetrics$clade == 2, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)
  
}



plot(c(0.8,2.6), range(c(scaledMetrics3$clade1, scaledMetrics3$clade2), na.rm = T), type = 'n', xaxt = 'n', xlab = "", ylab = "Scaled tree metric", las = 1, cex.lab = 1.3, main = "", cex.main = 1.3)
axis(1, at = 1:2, labels = c('Swifts', 'Hummingbirds'), tck = -.01, cex.axis = 1.3)

for (m in scaledMetrics3$metric) {
  
  points(x = c(1, 2), y = scaledMetrics3[scaledMetrics3$metric == m, 2:3], 
         type = 'b', col = scaledMetrics3$metricColor[scaledMetrics3$metric == m], lwd = 3)
  
  text(2.1, scaledMetrics3$clade2[scaledMetrics3$metric == m], m, adj = 0)
  
}
#text(.65, max(c(scaledMetrics3$clade1, scaledMetrics3$clade2), na.rm = T), "Metrics", cex = 1.5)
mtext("C) Tree metrics correlated with env. filtering", 3, adj = 0, cex = 1.3, bty='n', line = 1)


# Inferred environmental filtering parameters for different models
envParams = filter(scaledPredictions2, pair == 4, predictor == 'env1')

plot(c(0.8,2.6), range(envParams$scaledValue), type = 'n', xaxt = 'n', xlab = "", ylab = "low <-- Inferred env filtering --> high", las = 1, cex.lab = 1.3, main = "", cex.main = 1.3)
axis(1, at = 1:2, labels = c('Swifts', 'Hummingbirds'), tck = -.01, cex.axis = 1.3)

for (m in unique(envParams$model)) {
  
  points(envParams$clade[envParams$model==m], envParams$scaledValue[envParams$model==m], 
         type = 'b', col = envParams$color[envParams$model==m], lwd = 3)
  
  text(2.1, envParams$scaledValue[envParams$model==m & envParams$clade == 2], m, adj = 0)
  
}
#text(.65, max(envParams$scaledValue, na.rm = T), "Models", cex = 1.5)
mtext("D) Inference from simulation models", 3, adj = 0, cex = 1.3, bty='n', line = 1)

dev.off()




########################################################################################
# Example 3: Lagomorpha vs Rodentia and niche conservatism

clade1 = read.tree('trees/empirical/pair5_lagomorpha sister clade rodentia.txt')
clade2 = read.tree('trees/empirical/pair5_lagomorpha.txt')

pdf('figures/lagomorpha_nic_non-S-correlated_metrics.pdf', height = 8, width = 10)
par(mfcol = c(2, 2), mar = c(3, 4, 3, 1), oma = c(1, 0, 3, 0))
plot(clade1, show.tip.label = F, main = "", cex.main = 1.5)
mtext("A) Rodentia", 3, adj = 0, cex = 1.2, bty='n', line = 1)
plot(clade2, show.tip.label = F, main = "", cex.main = 1.5)
mtext("B) Lagomorpha", 3, adj = 0, cex = 1.2, bty='n', line = 1)

# Metrics associated with niche conservatism, at least 80% agreement (excluding those correlated with S)

nicMetrics = c('MGL_peakedness')


# Get scaled values of tree metrics relative to the distribution of metrics from all simulated trees
# In this case it seemed easier/more straightforward to calculate quantiles
# (in part because the distribution for many of these simulated metrics is highly skewed)
metColors = brewer.pal(length(nicMetrics), 'Blues')

scaledMetrics4 = data.frame(metric = nicMetrics, 
                            clade1 = rep(NA, length(nicMetrics)), 
                            clade2 = rep(NA, length(nicMetrics)),
                            metricColor = metColors[1:length(nicMetrics)])

for (d in nicMetrics) {
  
  scaledMetrics4$clade1[scaledMetrics4$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 5 & empiricalMetrics$clade == 1, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)
  
  scaledMetrics4$clade2[scaledMetrics4$metric == d] = sum(empiricalMetrics[empiricalMetrics$pair== 5 & empiricalMetrics$clade == 2, d] > processDFmetrics[, d], na.rm = T)/nrow(processDFmetrics)
  
}



plot(c(0.8,2.6), range(c(scaledMetrics4$clade1, scaledMetrics4$clade2), na.rm = T), type = 'n', xaxt = 'n', xlab = "", ylab = "Scaled tree metric", las = 1, cex.lab = 1.3, main = "", cex.main = 1.3)
axis(1, at = 1:2, labels = c('Rodentia', 'Lagomorpha'), tck = -.01, cex.axis = 1.3)

for (m in scaledMetrics4$metric) {
  
  points(x = c(1, 2), y = scaledMetrics4[scaledMetrics4$metric == m, 2:3], 
         type = 'b', col = scaledMetrics4$metricColor[scaledMetrics4$metric == m], lwd = 3)
  
  text(2.1, scaledMetrics4$clade2[scaledMetrics4$metric == m], m, adj = 0)
  
}
#text(.65, max(c(scaledMetrics4$clade1, scaledMetrics4$clade2), na.rm = T), "Metrics", cex = 1.5)
mtext("C) Tree metrics correlated with niche conservatism", 3, adj = 0, cex = 1.2, bty='n', line = 1)


# Inferred niche conservatism parameters for different models
nicParams = filter(scaledPredictions2, pair == 5, predictor == 'nic1')

plot(c(0.8,2.6), range(nicParams$scaledValue), type = 'n', xaxt = 'n', xlab = "", ylab = "low <- Inferred niche conservatism ->  high", las = 1, cex.lab = 1.3, main = "", cex.main = 1.3)
axis(1, at = 1:2, labels = c('Rodentia', 'Lagomorpha'), tck = -.01, cex.axis = 1.3)

for (m in unique(nicParams$model)) {
  
  points(nicParams$clade[nicParams$model==m], nicParams$scaledValue[nicParams$model==m], 
         type = 'b', col = nicParams$color[nicParams$model==m], lwd = 3)
  
  text(2.1, nicParams$scaledValue[nicParams$model==m & nicParams$clade == 2], m, adj = 0)
  
}
#text(.65, max(nicParams$scaledValue, na.rm = T), "Models", cex = 1.5)
mtext("D) Inference from simulation models", 3, adj = 0, cex = 1.2, bty='n', line = 1)

dev.off()
