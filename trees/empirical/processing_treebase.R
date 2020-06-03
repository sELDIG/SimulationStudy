# Modified from Dryad files for treebase
# https://datadryad.org/bitstream/handle/10255/dryad.96064/data-set_empirical_trees.zip?sequence=1

library(TreeSim)
library(apTreeshape)
load("trees/empirical/treebase_input.rda")

numbtrees <- length(treebase)
treeswithlength <- NULL # creates a list of trees that can be used for length statistics
####
####
#########
####
####
treeOutput = data.frame(model = NA, simID = NA, S = NA, tree.length = NA, PD = NA, gamma = NA, 
                        beta = NA, Colless = NA, Sackin = NA, Yule.PDA.ratio = NA, MRD = NA, 
                        VRD = NA, PSV = NA, mean.Iprime = NA, MPD = NA, VPD = NA, 
                        MGL_principal_eigenvalue = NA, MGL_asymmetry = NA,
                        MGL_peakedness = NA, MGL_eigengap = NA, nLTT_stat = NA)

write.table(treeOutput, 'treebase_empirical_output.txt', sep = '\t', row.names = FALSE)


for (i in 1:numbtrees)
{
	print(i)	
	{
	if ((class(treebase[[i]]) == "phylo")  )
	{
		
		{
		if (treebase[[i]]$Nnode >= 20 )	
		{
			#{
			#if (is.ultrametric(treebase[[i]]) == TRUE)
			#{
			  tree = treebase[[i]]
			  
			  model = 'treebase'
			  simID = i
			  
			  metrics = tryCatch({
			    treeMetrics(tree)
			  }, error = function(e) {
			    metrics =  data.frame(model = 'treebase', simID = NA, S = NA, tree.length = NA, PD = NA, gamma = NA, 
			                          beta = NA, Colless = NA, Sackin = NA, Yule.PDA.ratio = NA, MRD = NA, 
			                          VRD = NA, PSV = NA, mean.Iprime = NA, MPD = NA, VPD = NA, 
			                          MGL_principal_eigenvalue = NA, MGL_asymmetry = NA,
			                          MGL_peakedness = NA, MGL_eigengap = NA, nLTT_stat = NA)
			  })
			  
			  sink('treebase_empirical_output.txt', append = TRUE)
			  cat(paste(model, 
			            simID,
			            metrics$S, 
			            metrics$tree.length,
			            metrics$PD, 
			            metrics$gamma,
			            metrics$beta, 
			            metrics$Colless, 
			            metrics$Sackin,
			            metrics$Yule.PDA.ratio, 
			            metrics$MRD, 
			            metrics$VRD, 
			            metrics$PSV, 
			            metrics$mean.Iprime, 
			            metrics$MPD, 
			            metrics$VPD,
			            metrics$MGL_principal_eigenvalue, 
			            metrics$MGL_asymmetry, 
			            metrics$MGL_peakedness, 
			            metrics$MGL_eigengap, 
			            metrics$nLTT_stat,
			            '\n',
			            sep = '\t'))
			  sink()
			  
			}
			}
		}
		}
	}
	}
}
			  
