
# Package install and info

# library(Rcpp)
# library(devtools)
# install_url("https://dl.dropboxusercontent.com/s/wub79etldv8fry2/phylosim_0.3.1.tar.gz")


# ?PhyloSim
# browseVignettes("PhyloSim")

library(PhyloSim)
library(ape)

# And here the code to create new phylogenies. See options of the function to turn various things on and off


parameterList <- expand.grid(dispersal = c(0,1), density = c(0,1), environment = c(0,1), speciationRate = c(2,5), fission = c(1,2,3), protracted = c(0,2), replicate = c(1,2,3,4,5))

parameterList$id = 1:nrow(parameterList)

write.csv(parameterList, file = "fh-parameters.csv")

simulations = list()

for(i in 1:nrow(parameterList)){
  par <- createCompletePar(x = 100, y = 100, dispersal = parameterList$dispersal[i] , runs = c(5000), density = parameterList$density[i], environment = parameterList$environment[i], specRate = parameterList$speciationRate[i], protracted = parameterList$protracted[i]  )

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

