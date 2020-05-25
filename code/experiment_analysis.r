# Join parameter files indicating high, medium, and low treatment levels to tree output.

library(dplyr)

metrics = read.table('output.txt', header = T, sep = '\t', stringsAsFactors = FALSE)

models = unique(metrics$model)

joinedOutput = data.frame(model = NA, simID = NA, S = NA, log10S = NA, tree.length = NA, PD = NA, Bamma = NA, 
                          Beta = NA, Colless = NA, Sackin = NA, Yule.PDA.ratio = NA, MRD = NA, 
                          VRD = NA, PSV = NA, mean.Iprime = NA, MPD = NA, VPD = NA, 
                          MGL_principal_eigenvalue = NA, MGL_asymmetry = NA,
                          MGL_peakedness = NA, MGL_eigengap = NA, nLTT_stat = NA)

for (m in models) {
  param = read.csv(paste('parameters/', m, '_parameters.csv', sep = ''), header = T, stringsAsFactors = FALSE)
  
  treatments = param[, names(param) %in% c('model', 'simID', 'env', 'nic', 'dis', 'mut', 'tim')]
  
  # Add experiment columns if missing
  if(! 'env' %in% names(treatments)) treatments$env = NA
  if(! 'nic' %in% names(treatments)) treatments$nic = NA
  if(! 'dis' %in% names(treatments)) treatments$dis = NA
  if(! 'mut' %in% names(treatments)) treatments$mut = NA
  if(! 'tim' %in% names(treatments)) treatments$tim = NA
  
  # Order columns for consistency
  treatments = select(model, simID, env, nic, dis, mut, tim)
  
  # Join treatment columns to metrics output for a given model
  metricsSubset = filter(metrics, model == m) %>%
    left_join(treatments, by = c('model', 'simID'))
  
  joinedOutput = rbind(joinedOutput, metricsSubset)
  
}
joinedOutput = joinedOutput[-1, ]

