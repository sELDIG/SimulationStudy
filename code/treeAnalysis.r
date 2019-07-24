# Script for calculating metrics on all simulated phylogenies and conducting PCA
library(ape)
library(stringr)
library(dplyr)
library(geiger)

source('code/treeMetrics.R')

treeOutput = metricsForManyTrees()


treeMetricsPCA = function(treeOutput) {
  
}
pca = treeOutput %>%
  select(S, gamma.stat, Colless:mean.Iprime) %>%
  princomp()

# Function that joins 
classifyAcrossModels = function(treeOutput, crossModelSimTable) {
  classified = left_join(treeOutput, crossModelSimTable, by = c('model', 'simID'))
  return(classified)
}

classifyWithinModel = function(treeOutput, withinModelParametersTable) {
  classified = left_join(treeOutput, withinModelParametersTable, by = 'simID')
  return(classified)
}


pcaPlot = function(pca,                 # dataframe with model, simID, and PC scores
                   xscore = 1,          # PC score plotted on the x-axis
                   yscore = 2,          # PC score plotted on the y-axis
                   colorBy 'turquoise', # any variable/parameter name by which to color points
                   pchBy = 16           # categorical variable/parameter name by which to specify point shape 
                   ) {
  
  if (class(treeOutput) == 'numeric') {
    
    shades <- rainbow(130)[100:1]
    
    percents <- as.integer(cut(var, 100, 
                               include.lowest = TRUE, ordered = TRUE))
    fills <- shades[percents]
    
  }
  par(mar = c(4, 4, 1, 1))
  plot(pca$scores[,xscore], pca$scores[,yscore], col = as.character(treeOutput[,colorBy]), cex = 2, 
       pch = treeOutput$pch, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-max(abs(range(pca$scores[,yscore]))), max(abs(range(pca$scores[,yscore])))),
       xlim = c(-max(abs(range(pca$scores[,xscore]))), max(abs(range(pca$scores[,xscore])))))
  legend("topright", 
         legend = c("Pontarp", "DAISIE", "HS no zero-sum", "HS energy grad", "HS specn grad", "HS disturb grad"),
         col = c("orangered", "darkblue", rep("turquoise", 4)), pch = c(17, 17, 16, 17, 15, 1))
  
  par(new = TRUE)
  plot(1, 1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       xlim = c(-1, 1), ylim = c(-1, 1))
  text(pca$loadings[, xscore], pca$loadings[, yscore], row.names(pca$loadings))  
}
