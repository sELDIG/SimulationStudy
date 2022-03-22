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

