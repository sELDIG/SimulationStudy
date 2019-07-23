# Script for calculating metrics on all simulated phylogenies and conducting PCA
library(ape)
library(stringr)
library(dplyr)
library(geiger)

source('code/treeMetrics.R')

metricsForManyTrees = function(treefiles = NULL, treeOutput = NULL, minimumTreeSize = 20,
                               write = FALSE, fileOut) {

  if(is.null(treefiles)) {
    treefiles = list.files('trees')[grepl(".tre", list.files("trees"))]
  }  
  if(is.null(treeOutput)) {
    treeOutput = data.frame(model = NA, simID = NA, S = NA, gamma.stat = NA, beta.stat = NA, Colless = NA, 
                            Sackin = NA, shape.stat = NA, MRD = NA, VRD = NA, PSV = NA, mean.Iprime = NA)
  }
  
  for (treefile in treefiles) {
    
    treeIn = read.tree(paste("trees/", treefile, sep = ""))
    tree = drop.fossil(treeIn)
    
    if(tree$Nnode + 1 >= minimumTreeSize) {
      model = str_extract(treefile, "^[A-Za-z]*")
      simID = str_extract(treefile, "[0-9]+")
      
      print(treefile)
      metrics = treeMetrics(tree)
      
      treeOutput = rbind(treeOutput,
                         data.frame(model = model, simID = simID, S = metrics$S, gamma.stat = metrics$gamma.stat,
                                    beta.stat = metrics$beta.stat, Colless = metrics$Colless, Sackin = metrics$Sackin,
                                    shape.stat = metrics$shape.stat, MRD = metrics$MRD, VRD = metrics$VRD, 
                                    PSV = metrics$PSV, mean.Iprime = metrics$mean.Iprime,
                                    MGL_principal_eigenvalue = metrics$MGL_principal_eigenvalue, 
                                    MGL_asymmetry = metrics$MGL_asymmetry, 
                                    MGL_peakedness = metrics$MGL_peakedness, MGL_eigengap = metrics$MGL_eigengap))
    } else {
      print(paste(treefile, "skipped -- not enough species"))
    }
    
    if(write) {
      write.csv(treeOutput, paste("treeOutput_", fileOut, "_", Sys.Date(), ".csv", sep = ""), row.names = F)
    }
    
  }  
  treeOutput = filter(treeOutput, !is.na(model), !is.na(simID))
  return(treeOutput)    
}

treeOutput = metricsForManyTrees()

pca = treeOutput %>%
  select(S, gamma.stat, Colless:mean.Iprime) %>%
  princomp()


pcaPlot = function(xscore = 1, yscore = 2) {
  par(mar = c(4, 4, 1, 1))
  plot(pca$scores[,xscore], pca$scores[,yscore], col = as.character(treeOutput$col), cex = 2, 
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
