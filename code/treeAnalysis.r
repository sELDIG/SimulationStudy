library(ape)
library(stringr)
library(dplyr)
library(geiger)
library(lessR)
library(gsheet)
library(dplyr)

source('code/treeMetrics.R')
source('code/pcaFunctions.r')

### Analysis steps

# 0 - read in trees and calculate metrics using metricsForManyTrees()
# metricsForManyTrees(treefiles, fileOut = "treeOutput.txt", append = F)

# 1 - read in output
treeOutput = read.table("treeOutput.txt", header = T, sep = '\t', stringsAsFactors = FALSE)

# 2 - decide which vars to include in PCA
varsForPCA = c("PD", "Gamma", "Beta", "Colless", "Sackin", "Yule.PDA.ratio", "MRD",
               "VRD", "PSV", "mean.Iprime", "MPD", "MGL_principal_eigenvalue",
               "MGL_asymmetry", "MGL_peakedness", "MGL_eigengap", "nLTT_stat")

varsIndependent = c("Gamma", "Beta", "MRD", "MPD", "MGL_peakedness")


# 3 - run PCA (default on all models together)
# If you want to run PCA on just one or a few models, specify their abbreviation(s) with 'models'
# e.g., models = 'hs'
pcaOutput = treeMetricsPCA(treeOutput, models = 'all', vars = varsForPCA)

# 4 - Visualize between simulation models; 
#     choose which model classification variables to visualize by color or symbol, 
#     and on which PC axes
betweenModelPCAPlot(pcaOutput, xscore = 1, yscore = 2, colorBy = "model", pchBy = "ModelFamily")

# 5 - Visualize within simulation models;
#     choose which model classification variables to visualize by color or symbol
withinModelPCAPlot(pcaOutput, "ra", xscore = 2, yscore = 3, colorBy = 'Dispersal', pchBy = 'Founder')



# 6 - Empirical data; all clades/subclades of mammals and amphibians greater than 20 species
mammalfiles = list.files('trees/empirical/mammal_clades/')

mammalOut = metricsForManyTrees(mammalfiles, fileOut = "mammalMetrics.txt", append = FALSE)
mammalOutput = read.table("mammalMetrics.txt", sep = '\t', header = T, fill = TRUE) %>%
  filter(!is.na(S), S > 20)

