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
# 1 - read in output
treeOutput = read.table("treeOutput.txt", header = T, sep = '\t')

# 2 - decide which vars to include in PCA
varsForPCA = c("PD", "gamma", "beta", "Colless", "Sackin", "Yule.PDA.ratio", "MRD",
               "VRD", "PSV", "mean.Iprime", "MPD", "MGL_principal_eigenvalue",
               "MGL_asymmetry", "MGL_peakedness", "MGL_eigengap", "nLTT_stat")

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





