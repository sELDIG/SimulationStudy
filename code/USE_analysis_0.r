# Analyzing Uniform Sampling Experiment trees

# Load libraries
library(gsheet)
library(dplyr)
library(stringr)
library(ape)
library(geiger)
library(lessR)

source('code/pcaFunctions.r')
source('code/experiment_analysis_functions.r')
source('code/treeMetrics.r')

# Parameter key across simulations; 
# url <- 'https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=1171496897'
# paramKey <- gsheet2tbl(url)
# write.csv(paramKey, 'experiments/uniform_sampling_experiment/simulation_parameters_key.csv', row.names = F)

paramKey = read.csv('experiments/uniform_sampling_experiment/simulation_parameters_key.csv')

# Align parameter values with the relevant experiment, create a dataframe with columns for
# model, model2 (with scenario suffix if appropriate), simID, and columns for the parameter
# names and parameter values associated with 5 processes: env, dis, nic, mut, com

modelList = c('ca', 'fh', 'gen', 'hs', 'pontarp', 've', 'xe')

for (m in modelList) {
  
  if (!exists("processDF")) {
    processDF = alignParametersWithProcesses(m)
  } else {
    processDF = rbind(processDF, alignParametersWithProcesses(m))
  }
  
}



# Join tree metrics to the aligned parameter-process dataframe: for analysis, use **processDFmetrics**
# Script for calculating tree metrics for trees that have not already been analyzed

# Read in existing output file
treeOutput = read.table('USE_treeOutput.txt', sep = '\t', header = T)

analyzedTrees = paste(treeOutput$model, "_", treeOutput$simID, ".tre", sep = "")

# Get list of tree files in USE experiment folder
allTreeFiles = list.files('trees/uniform_sampling_experiment')[grepl(".tre", list.files('trees/uniform_sampling_experiment'))]

treesToRun = allTreeFiles[!allTreeFiles %in% analyzedTrees]

# Calculate tree metrics for remaining trees and append
metricsForManyTrees(treefiles = treesToRun, minimumTreeSize = 5, fileOut = 'USE_treeOutput.txt', append = TRUE, 
                    treedir = 'trees/uniform_sampling_experiment')


# Join process-parameter linkage dataframe to tree output
processDFmetrics = left_join(processDF, treeOutput, by = c('model', 'simID'))

write.csv(processDFmetrics, 'experiments/uniform_sampling_experiment/process_parameter_values_and_tree_metrics.csv', 
          row.names = F)

