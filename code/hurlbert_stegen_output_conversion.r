# Convert original raw simulation output from Hurlbert & Stegen model
# to 

library(ape)
library(dplyr)
library(geiger)

sim_vector = c()

for (i in sim_vector) {
  
  tree = read.tree(paste('raw_sim_output/sim', i, '_out/SENC_phylo_sim', i, '.tre', sep = ''))
  #tree = drop.extinct(tree)
  write.tree(tree, paste('z:/git/simulationstudy/trees/hs_', i, '.tre', sep = ''))
  
  all.pops = read.csv(paste('raw_sim_output/sim', i, '_out/SENC_all.pops_sim', i, '.csv', sep = ''))
  
  # Species distributions at the end of the simulation
  extant.pops = all.pops %>%
    filter(extant == 1)
  
  sites = extant.pops %>%
    select(spp.name, region) %>%
    distinct()
  write.csv(sites, paste('z:/git/simulationstudy/sitexspecies/hs_', i, '_sites.csv', sep = ''))
  
  traits = extant.pops %>%
    select(spp.name, env.opt) %>%
    distinct()
  write.csv(sites, paste('z:/git/simulationstudy/traits/hs_', i, '_traits.csv', sep = ''))
  
}