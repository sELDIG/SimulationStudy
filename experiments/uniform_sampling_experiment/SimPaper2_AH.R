<<<<<<< HEAD
=======
# This script 
# 1) examines simulated trees and metrics and relates them to simulation parameter values.
# 2) identifies which process-parameters are most correlated with which tree metrics
# 3) conducts random forest models for describing which metrics are most predictive of parameter values
# 4) inverts the random forest models to infer estimated parameter values given tree shape metrics for several pairs of sister clades

library(RColorBrewer)
library(dplyr)

# Read in simulated tree data
>>>>>>> 1614d5e1bc64cd80d94f900ffebc7f681457894d
simuData <- read.csv("./experiments/uniform_sampling_experiment/process_parameter_values_and_tree_metrics_sign_corrected.csv", stringsAsFactors = T)
simuDataKey <- read.csv("./experiments/uniform_sampling_experiment/simulation_parameters_key.csv")


head(simuDataKey)

statistics = colnames(simuData)[25:43]
statisticsIndices = 25:43

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

 

<<<<<<< HEAD
image.real <- function(mat, xCol = c("blue", "white", "white", "red"), range = c(-1,1)) { 
  mat <- t(mat)[,nrow(mat):1]
  fields::image.plot(mat, axes = FALSE, zlim = range, 
                     col = colorRampPalette(xCol)(30))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat))
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat), las = 2)
=======
image.real <- function(mat, xCol = c("blue", "white", "white", "red"), 
                       range = c(-1,1), x.labels = rownames(mat), y.labels = colnames(mat)) { 
  mat <- t(mat)[,nrow(mat):1]
  fields::image.plot(mat, axes = FALSE, zlim = range, 
                     col = colorRampPalette(xCol)(30))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = x.labels)
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = y.labels, las = 2)
>>>>>>> 1614d5e1bc64cd80d94f900ffebc7f681457894d
  box() 
}


# 2-panel figure
par(mfrow = c(1, 2), mar = c(3,11,3,3))

# Values: % agreement, Color: average sign
<<<<<<< HEAD
image.real(signMat) 
=======
image.real(signMat, x.labels = c('env', 'dis', 'nic', 'mut', 'com')) 
>>>>>>> 1614d5e1bc64cd80d94f900ffebc7f681457894d
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
<<<<<<< HEAD
image.real(signMat) 
=======
image.real(signMat, x.labels = c('env', 'dis', 'nic', 'mut', 'com')) 
>>>>>>> 1614d5e1bc64cd80d94f900ffebc7f681457894d
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


<<<<<<< HEAD

=======
########################################################################################
# Random Forests for predicting parameter values from tree metrics
>>>>>>> 1614d5e1bc64cd80d94f900ffebc7f681457894d

library(ranger)

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

<<<<<<< HEAD
=======
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


# Random forests for each model, saving predictions of parameter values from random forests for each empirical tree
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
         clade = rep(c(1, 2), 174),
         scaledValue = apply(scaledPredictions[, c('value', 'minVal', 'maxVal')], 1, 
                             function(x) scaleFunction(value = x[1], min = x[2], max = x[3]))) %>%
  left_join(modelColors, by = 'model')




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

dispMetrics = c('VRD', 'mean.Iprime', 'Colless', 'Sackin', 'Yule.PDA.ratio', 'MRD', 'MGL_principal_eigenvalue')


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

>>>>>>> 1614d5e1bc64cd80d94f900ffebc7f681457894d
