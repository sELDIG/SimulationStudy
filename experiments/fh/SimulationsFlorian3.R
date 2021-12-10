# Started 2.6.21
# New simulations 

library(PhyloSim)
library(ape)

set.seed(123)

n = 2000

dispersalGlobal = sample(c(T,F), n, replace = T, prob = c(0.2, 0.8))

parameterList <- data.frame(
  model = rep("fh", n),
  simID = 1:n,
  dispersal = ifelse(dispersalGlobal, 0, runif(n, 1, 5)),
  density = runif(n, 0, 1),
  environment = runif(n, 0, 1),
  speciationRate = runif(n, 0.5, 5),
  fission = sample.int(3, n, replace = T)
)

parameterList$scenario = parameterList$fission
parameterList$dispersalGlobal = dispersalGlobal

expName = "experiment6"
dir.create(expName)
dir.create(paste(expName, "/plots", sep = ""))
dir.create(paste(expName, "/trees", sep = ""))
write.csv(parameterList, file = paste(expName, "/fh-parameters.csv", sep = ""))
write.csv(parameterList, file = "../../trees/uniform_sampling_experiment/fh_USE_parameters.csv", row.names = F)


pars = list()

for(i in 1:n){
  pars[[i]] <- createCompletePar(x = 100, y = 100, 
                                 dispersal = parameterList$dispersal[i] , 
                                 runs = c(50000), 
                                 density = parameterList$density[i], 
                                 environment = parameterList$environment[i], 
                                 specRate = parameterList$speciationRate[i])
}


simulations <- runSimulationBatch(pars, parallel = 10)

save(simulations, file = paste(expName, "/simulations.Rdata", sep = ""))

for(i in 1:n){
  
  pdf(file = paste(expName, "/plots/", i, ".pdf", sep=""))
  
  plotSpatialPhylo(simulations[[i]])
  dev.off()
  phylogeny <- simulations[[i]]$Output[[1]]$phylogeny
  
  extantPhylogeny <- drop.fossil(phylogeny)
  
  write.tree(extantPhylogeny, file =  paste(expName, "/trees/fh_", i, ".tre", sep=""))  
}










