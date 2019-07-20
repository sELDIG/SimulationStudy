# Script for calculating metrics on all simulated phylogenies and conducting PCA
library(ape)
library(stringr)
library(dplyr)

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
    
    tree = read.tree(paste("trees/", treefile, sep = ""))
    
    if(tree$Nnode + 1 >= minimumTreeSize) {
      model = str_extract(treefile, "^[A-Za-z]*")
      simID = str_extract(treefile, "[0-9]+")
      
      print(treefile)
      metrics = treeMetrics(tree)
      
      treeOutput = rbind(treeOutput,
                         data.frame(model = model, simID = simID, S = metrics$S, gamma.stat = metrics$gamma.stat,
                                    beta.stat = metrics$beta.stat, Colless = metrics$Colless, Sackin = metrics$Sackin,
                                    shape.stat = metrics$shape.stat, MRD = metrics$MRD, VRD = metrics$VRD, 
                                    PSV = metrics$PSV, mean.Iprime = metrics$mean.Iprime))
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



