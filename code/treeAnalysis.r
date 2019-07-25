# Script for calculating metrics on all simulated phylogenies and conducting PCA
library(ape)
library(stringr)
library(dplyr)
library(geiger)
library(lessR)
library(gsheet)
library(dplyr)

source('code/treeMetrics.R')

# Reads in model classification codes and assigns colors and pchs
url = "https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=2047946073"

modelClassification = gsheet2tbl(url)




# Function that conducts PCA on treeOutput set of specified tree metrics
# -- vars is a vector of column names to include, default is all vars
treeMetricsPCA = function(treeOutput, vars = 'all') {

  if (vars == 'all') {
    vars = names(treeOutput[, 3:ncol(treeOutput)])
  }
  outputSubsetNoNAs = na.omit(treeOutput[!names(treeOutput) == "VPD"])  
  pc = princomp(outputSubsetNoNAs[, names(outputSubsetNoNAs) %in% vars], cor = TRUE)
  pcaScores = cbind(outputSubsetNoNAs[, c("model", "simID")], pc$scores) 
  
  return(list(pcaScores = pcaScores, pcaLoadings = pc$loadings))
}


classifyWithinModel = function(treeOutput, withinModelParametersTable) {
  classified = left_join(treeOutput, withinModelParametersTable, by = 'simID')
  return(classified)
}


betweenModelPCAPlot = function(pcaOutput,          # dataframe with model, simID, and PC scores
                               xscore = 1,        # PC score plotted on the x-axis
                               yscore = 2,        # PC score plotted on the y-axis
                               colorBy = 'model', # any variable/parameter name by which to color points
                               pchBy = 16         # categorical variable/parameter name by which to specify point shape 
                               ) {
  require(gsheet)
  url = "https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=2047946073"
  
  modelClassification = gsheet2tbl(url)
  
  pcaScores = pcaOutput$pcaScores
  pcaLoadings = pcaOutput$pcaLoadings
  
  colors = c('turquoise', 'orangered', 'yellow2', 'darkblue', 'limegreen', 'magenta', 'blue', 'black')
  pch = c(15, 16, 17, 18, 1, 7, 8, 10)

  colorCode = data.frame(val = unique(modelClassification[, colorBy]), 
                         color = colors[1:nrow(unique(modelClassification[, colorBy]))])
  colorCode$color = as.character(colorCode$color)
  names(colorCode)[1] = colorBy
  
  pchCode = data.frame(val = unique(modelClassification[, pchBy]),
                       pch = pch[1:nrow(unique(modelClassification[, pchBy]))])
  names(pchCode)[1] = pchBy
  
  plotOutput = left_join(pcaScores, modelClassification, by = c("model", "simID")) %>%
    left_join(colorCode, by = unname(colorBy)) %>%
    left_join(pchCode, by = unname(pchBy))

  plot(plotOutput[,paste("Comp.", xscore, sep = "")], plotOutput[, paste("Comp.", yscore, sep = "")], 
       col = plotOutput$color, cex = 2, 
       pch = plotOutput$pch, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")]))), 
                max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")])))),
       xlim = c(-max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")]))), 
                max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")])))))
  legend("topleft", 
         legend = c(toupper(colorBy), colorCode[,1]), bty = "n",
         col = c('white', colorCode[,2]), pch = 16)
  
  legend("topright", 
         legend = c(toupper(pchBy), pchCode[,1]), bty = "n",
         col = c('white', rep('black', nrow(pchCode))), pch = c(16, pchCode[,2]))
  
  par(new = TRUE)
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       xlim = c(-1, 1), ylim = c(-1, 1))
  
  # Only print the top variable loadings since it otherwise becomes messy
  mainLoadings = pcaLoadings[order(pcaLoadings[, xscore], decreasing = TRUE) <= 4 | 
                               order(pcaLoadings[,yscore], decreasing = TRUE) <= 4, 
                             c(xscore, yscore)]
  text(mainLoadings[, 1], mainLoadings[, 2], row.names(mainLoadings))  
}


# Analysis
treeOutput = read.table("treeOutput.txt", header = T, sep = '\t')

varsForPCA = c("PD", "gamma", "beta", "Colless", "Sackin", "Yule.PDA.ratio", "MRD",
               "VRD", "PSV", "mean.Iprime", "MPD", "MGL_principal_eigenvalue",
               "MGL_asymmetry", "MGL_peakedness", "MGL_eigengap", "nLTT_stat")

pcaOutput = treeMetricsPCA(treeOutput, vars = varsForPCA)

par(mfrow = c(1,1))
pcaPlot(pca, xscore = 1, yscore = 2)
pcaPlot(pca, xscore = 1, yscore = 3)
pcaPlot(pca, xscore = 2, yscore = 3)
pcaPlot(pca, xscore = 2, yscore = 4)

varCor = cor(treeOutput[,!names(treeOutput) %in% c('model', 'simID', 'VPD')], use = "pairwise.complete.obs")
varCor2 = corReorder(varCor, bottom = 6, right = 6, diagonal_new = FALSE)



#if (class(treeOutput) == 'numeric') {

# shades <- rainbow(130)[100:1]

#percents <- as.integer(cut(var, 100, 
#include.lowest = TRUE, ordered = TRUE))
#fills <- shades[percents]

#}

