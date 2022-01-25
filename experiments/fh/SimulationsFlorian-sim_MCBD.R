# Started 2.6.21
# New simulations 

library(ape)
library(RPANDA)
library(parallel)

set.seed(123)

n = 500


parameterList <- data.frame(
  model = rep("la", n), # Leandro Aristide
  simID = 1:n,
  density = runif(n, 0.01, 0.1), # pars[8] corresponds to alpha1, the competition effect on extinction (competition strength)
  speciationRate = runif(n, 0.1, 0.5) # pars[1] corresponds to lambda1, the speciation intitation rate
)


expName = "experiment-sim_MCBD"
dir.create(expName)
dir.create(paste(expName, "/plots", sep = ""))
dir.create(paste(expName, "/trees", sep = ""))
# choosing here la for Leandro Aristide 
write.csv(parameterList, file = paste(expName, "/la-parameters.csv", sep = ""))
write.csv(parameterList, file = "../../trees/uniform_sampling_experiment/la_USE_parameters.csv", row.names = F)


pars = list()

for(i in 1:n){
  lambda1 = parameterList$speciationRate[i]
  tau0 = 0.01
  beta = 0.6
  mu0 = 0.5
  mubg = 0.01
  mui0 = 0.8
  muibg = 0.02
  alpha1 = alpha2 = parameterList$density[i]
  sig2 = 0.5
  m = 20
  
  pars[[i]] <- c(lambda1, tau0, beta, mu0, mubg,mui0, muibg, alpha1, alpha2, sig2, m)
}


getResults <- function(par){
  out = NA
  #out = RPANDA::sim_MCBD(par, age.max=100, step.size=0.1, plot = F)
  out = try(RPANDA::sim_MCBD(par, age.max=100, step.size=0.1, plot = F), silent = TRUE)
  
  out = try(out$gsp_extant$tree, silent = T)
  return(out)
}


# simulations <- lapply(pars, getResults)

clust <- makeCluster(10)
simulations <- parLapplyLB(clust, pars, getResults)
stopCluster(clust)


save(simulations, file = paste(expName, "/simulations.Rdata", sep = ""))

for(i in 1:n){
  
  if(class(simulations[[i]]) == "phylo"){
    pdf(file = paste(expName, "/plots/", i, ".pdf", sep=""))
    
    plot(simulations[[i]])
    dev.off()
    write.tree(simulations[[i]], file =  paste(expName, "/trees/la_", i, ".tre", sep=""))      
  }
}



