simuData <- read.csv("./experiments/uniform_sampling_experiment/process_parameter_values_and_tree_metrics.csv", stringsAsFactors = T)
simuDataKey <- read.csv("./experiments/uniform_sampling_experiment/simulation_parameters_key.csv")



statistics = colnames(simuData)[-(1:24)]
statisticsIndices = 25:ncol(simuData)

predictorsIndices = c(14,16,18,20,22)
predictors = colnames(simuData)[c(14,16,18,20,22)]

signMat = matrix(nrow = length(statistics), ncol = length(predictors))
rownames(signMat) = statistics
colnames(signMat) = predictors
R2Mat = disAg = signMat

models = levels(simuData$model)


for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    R2 = sig = rep(NA, length(models))
    for(k in 1:length(models)){
      tmp = simuData[simuData$model == models[k],]
      x = scale(tmp[,statisticsIndices[i]])
      y = tmp[,predictorsIndices[j]]
      if(! all(is.na(y))){
      fit <- summary(lm(y ~ x))
      R2[k] = fit$r.squared
      sig[k] = ifelse(fit$coefficients[2,4] <0.05, sign(fit$coefficients[2,1]), 0)
      } else{
        R2[k] = NA
        sig[k] = NA
      }
    }
    R2Mat[i,j] = mean(R2, na.rm = T)
    signMat[i,j] = sum(sig, na.rm = T) / sum(!is.na(sig))
    disAg[i,j] = abs(sum(sig == 1, na.rm = T) - sum(sig == -1, na.rm = T)) / abs(sum(sig %in% c(1,-1), na.rm = T))
  }
}
 

image.real <- function(mat, xCol = c("blue", "white", "white", "red"), range = c(-1,1)) { 
  mat <- t(mat)[,nrow(mat):1]
  fields::image.plot(mat, axes = FALSE, zlim = range, 
                     col = colorRampPalette(xCol)(30))
  axis(1, at = seq(0, 1, length = nrow(mat)), labels = rownames(mat))
  axis(2, at = seq(0, 1, length = ncol(mat)), labels = colnames(mat), las = 2)
  box() 
}

par(mar = c(3,12,3,3))

image.real(signMat) 
for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices)-1), 
         labels = round(disAg[i,j], digits = 2))
  }
}

title(main = "Model agreement regarding effect direction")


image.real(R2Mat, range = c(0,1), xCol = c("beige", "firebrick2")) 
title(main = "R2 between parameter and tree metric")

for(i in 1:length(statisticsIndices)){
  for(j in 1:length(predictorsIndices)){
    text((j-1)/(length(predictorsIndices)-1),
         1-(i-1)/(length(statisticsIndices)-1), 
         labels = round(R2Mat[i,j], digits = 2))
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

