# Script for calculating metrics on all simulated phylogenies and conducting PCA
library(ape)
library(stringr)
library(dplyr)
library(geiger)
library(lessR)

source('code/treeMetrics.R')

varsIndependentOfS = c("PD", "gamma", "beta", "Colless", "Sackin", "Yule.PDA.ratio", "MRD",
                       "VRD", "PSV", "mean.Iprime", "MPD", "MGL_principal_eigenvalue",
                       "MGL_asymmetry", "MGL_peakedness", "MGL_eigengap", "nLTT_stat")

modelColorAssignment = function() {
  colors = c('turquoise', 'magenta', 'red', 'darkblue', 'limegreen', 'yellow2', 'blue', 'black')
  
  crossModelSimTable = data.frame(model = unique(treeOutput$model), 
                                  color = colors[1:length(unique(treeOutput$model))])
  crossModelSimTable$color = as.character(crossModelSimTable$color)
  
  return(crossModelSimTable)
}

treeMetricsPCA = function(treeOutput, vars = NULL) {

  crossModelSimTable = modelColorAssignment()
    
  if (is.null(vars)) {
    vars = names(treeOutput[, 3:ncol(treeOutput)])
  }
  outputSubsetNoNAs = na.omit(treeOutput[!names(treeOutput) == "VPD"])  
  pc = princomp(outputSubsetNoNAs[, names(outputSubsetNoNAs) %in% vars], cor = TRUE)
  pcaOutput = cbind(outputSubsetNoNAs[, c("model", "simID")], pc$scores) %>%
    left_join(crossModelSimTable, by = 'model')
  
  return(list(pca = pcaOutput, loadings = pc$loadings))
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


pcaPlot = function(pc,                 # dataframe with model, simID, and PC scores
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
  
  crossModelSimTable = modelColorAssignment()
  
  if ("color" %in% names(pc$pca)) { colorBy = pc$pca$color }
  if ("pch" %in% names(pc$pca)) { colorBy = pc$pca$pch }
  
  plot(pc$pca[,paste("Comp.", xscore, sep = "")], pc$pca[, paste("Comp.", yscore, sep = "")], 
       col = colorBy, cex = 2, 
       pch = pchBy, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-max(abs(range(pc$pca[,paste("Comp.", yscore, sep = "")]))), 
                max(abs(range(pc$pca[,paste("Comp.", yscore, sep = "")])))),
       xlim = c(-max(abs(range(pc$pca[,paste("Comp.", xscore, sep = "")]))), 
                max(abs(range(pc$pca[,paste("Comp.", xscore, sep = "")])))))
  legend("topleft", 
         legend = crossModelSimTable$model,
         col = crossModelSimTable$color, pch = pchBy)
  
  par(new = TRUE)
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       xlim = c(-1, 1), ylim = c(-1, 1))
  
  # Only print the top variable loadings since it otherwise becomes messy
  mainLoadings = pc$loadings[order(pc$loadings[, xscore], decreasing = TRUE) <= 4 | 
                               order(pc$loadings[,yscore], decreasing = TRUE) <= 4, 
                             c(xscore, yscore)]
  text(mainLoadings[, 1], mainLoadings[, 2], row.names(mainLoadings))  
}


# Analysis
treeOutput = read.table("treeOutput.txt", header = T, sep = '\t')

pca = treeMetricsPCA(treeOutput, vars = varsIndependentOfS)

par(mfrow = c(1,1))
pcaPlot(pca, xscore = 1, yscore = 2)
pcaPlot(pca, xscore = 1, yscore = 3)
pcaPlot(pca, xscore = 2, yscore = 3)
pcaPlot(pca, xscore = 2, yscore = 4)

varCor = cor(treeOutput[,!names(treeOutput) %in% c('model', 'simID', 'VPD')], use = "pairwise.complete.obs")
varCor2 = corReorder(varCor, bottom = 6, right = 6, diagonal_new = FALSE)
