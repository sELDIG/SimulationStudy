library(PhyloSim)
library(ape)

# And here the code to create new phylogenies. See options of the function to turn various things on and off

set.seed(123)

n = 5000

dispersalGlobal = sample(c(T,F), n, replace = T, prob = c(0.2, 0.8))

parameterList <- data.frame(
  dispersal = ifelse(dispersalGlobal, 0, runif(n, 1, 5)),
  density = runif(n, 0, 1),
  environment = runif(n, 0, 1),
  speciationRate = runif(n, 0.5, 5),
  fission = sample.int(3, n, replace = T),
  id = 1:n
)

parameterList$dispersalGlobal = dispersalGlobal

  

write.csv(parameterList, file = "fh-parameters-2.csv")

simulations = list()

for(i in 1:nrow(parameterList)){
  par <- createCompletePar(x = 100, y = 100, dispersal = parameterList$dispersal[i] , runs = c(5000), density = parameterList$density[i], environment = parameterList$environment[i], specRate = parameterList$speciationRate[i])

  new.simulation <- runSimulation(par)
  
  simulations[[i]] = new.simulation

  pdf(file = paste("./plots/", i, ".pdf", sep=""))
  
  plotSpatialPhylo(new.simulation)
  dev.off()
  phylogeny <- new.simulation$Output[[1]]$phylogeny
  
  extantPhylogeny <- drop.fossil(phylogeny)
  
  write.tree(extantPhylogeny, file =  paste("./trees/fh_", i, ".tre", sep=""))
}

save(simulations, file = "simulations.Rdata")


# the output is in simulation output. If several time steps are provided in runs, the output is a unnamed list for the time steps, and each entry will contain all state variables

# so, the phylogeny after 1000 time steps is in 



specRich(new.simulation)
sac(new.simulation)

