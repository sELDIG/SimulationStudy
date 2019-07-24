# Script for calculating metrics on all simulated phylogenies and conducting PCA
library(ape)
library(stringr)
library(dplyr)
library(geiger)

source('code/treeMetrics.R')

treeOutput = read.table("treeOutput.txt", header = T, sep = '\t')

colors = c('turquoise', 'magenta', 'orangered', 'darkblue', 'limegreen', 'yellow4', 'blue', 'black')

crossModelSimTable = data.frame(model = unique(treeOutput$model), 
                                color = colors[1:length(unique(treeOutput$model))])

varsIndependentOfS = c("PD", "gamma", "beta", "Colless", "Sackin", "Yule.PDA.ratio", "MRD",
                       "VRD", "PSV", "mean.Iprime", "MPD", "MGL_principal_eigenvalue",
                       "MGL_asymmetry", "MGL_peakedness", "MGL_eigengap", "nLTT_stat")

treeMetricsPCA = function(treeOutput, vars = NULL) {
  
  if (is.null(vars)) {
    vars = names(treeOutput[, 3:ncol(treeOutput)])
  }
  outputSubsetNoNAs = na.omit(treeOutput[!names(treeOutput) == "VPD"])  
  pc = princomp(outputSubsetNoNAs[, names(outputSubsetNoNAs) %in% vars], cor = TRUE)
  pcaOutput = cbind(outputSubsetNoNAs[, c("model", "simID")], pc$scores) %>%
    left_join(crossModelSimTable, by = 'model')
}

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
                   colorBy = 'turquoise', # any variable/parameter name by which to color points
                   pchBy = 16           # categorical variable/parameter name by which to specify point shape 
                   ) {
  
  #if (class(treeOutput) == 'numeric') {
    
   # shades <- rainbow(130)[100:1]
    
    #percents <- as.integer(cut(var, 100, 
                               #include.lowest = TRUE, ordered = TRUE))
    #fills <- shades[percents]
    
  #}
  par(mar = c(4, 4, 1, 1))
  plot(pca[,paste("Comp.", xscore, sep = "")], pca[, paste("Comp.", yscore, sep = "")], 
       col = as.character(pca$color), cex = 2, 
       pch = pchBy, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-max(abs(range(pca[,paste("Comp.", yscore, sep = "")]))), max(abs(range(pca[,paste("Comp.", yscore, sep = "")])))),
       xlim = c(-max(abs(range(pca[,paste("Comp.", xscore, sep = "")]))), max(abs(range(pca[,paste("Comp.", xscore, sep = "")])))))
  legend("topright", 
         legend = crossModelSimTable$model,
         col = crossModelSimTable$color, pch = pchBy)
  
  par(new = TRUE)
  plot(1, 1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       xlim = c(-1, 1), ylim = c(-1, 1))
  text(pca$loadings[, xscore], pca$loadings[, yscore], row.names(pca$loadings))  
}
