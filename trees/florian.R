# Hi guys, I wasn't really sure what to simulate, but we can do this any time. 
# Here are instructions to install the model as an R package

# install.packages(c("devtools","Rcpp")) # I case you don't have them installed

library(Rcpp)
library(devtools)
install_url("https://dl.dropboxusercontent.com/s/wub79etldv8fry2/phylosim_0.3.1.tar.gz")


?PhyloSim
browseVignettes("PhyloSim")

library(PhyloSim)

# And here the code to create new phylogenies. See options of the function to turn various things on and off

par <- createCompletePar(x = 100, y = 100, dispersal = "global" , 
runs = c(1000,10000),density = 1, environment =0.5, specRate = 1)
new.simulation <- runSimulation(par)

plotSpatialPhylo(new.simulation)


# the output is in simulation output. If several time steps are provided in runs, the output is a unnamed list for the time steps, and each entry will contain all state variables

# so, the phylogeny after 1000 time steps is in 

new.simulation$Output[[1]]$phylogeny

# and after 10.000 time steps, we have 

new.simulation$Output[[2]]$phylogeny
