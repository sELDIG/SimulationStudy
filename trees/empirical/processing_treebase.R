library(TreeSim)
library(apTreeshape)
###   load(treebase_input)

treebaselist <- list("age", "kind","ultrametrictrue","ultrametricfalse","problematicwithas.treeshape", "maxbt", "minbt", "Nnode","realleaves","realleavesforlength", "extinct", "samesp", "gammaStat", "colless", "sackin", "maxlik.betasplit")
numbtrees <- length(treebase)
treeswithlength <- NULL # creates a list of trees that can be used for length statistics
####
####
#########
####
####
i <- 1
# for mytree.sy.age
for (i in i:numbtrees)
{
	print(i)	
	{
	if ((class(treebase[[i]]) == "phylo")  )
	{
		
		{
		if (length(treebase[[i]]$edge.length) > 0 )	
		{
			{
			if (is.ultrametric(treebase[[i]]) == TRUE)
			{
				treebaselist$ultrametrictrue <- c(treebaselist$ultrametrictrue, i)
				treebaselist$age <- c(		treebaselist$maxbt , max(branching.times(treebase[[i]])))
				treebaselist$maxbt <- c(		treebaselist$maxbt , max(branching.times(treebase[[i]])))
				treebaselist$minbt <- c(		treebaselist$minbt, min(branching.times(treebase[[i]])))
				treebaselist$gammaStat <- c(		treebaselist$gammaStat, gammaStat(treebase[[i]]))
				treebaselist$realleavesforlength <- c(		treebaselist$realleavesforlength, length(treebase[[i]]$tip.label))
				
				treeswithlength[[length(treeswithlength)+ 1]] <- treebase[[i]]
					
			}
			else
			{
				treebaselist$ultrametricfalse <- c(treebaselist$ultrametricfalse, i)	
			}
			}
		}
		}
		{
		if((length(treebase[[i]]$tip.label) > 2)  & length(try(write.nexus(treebase[[i]] , file=paste("treebaseforbeta", numbtrees , ".nex", sep="")))) == 0 & length(try(as.treeshape(treebase[[i]]))) != 0 & treebase[[i]]$kind=="Species Tree")
		{
			{
			if  (length(as.treeshape(treebase[[i]])) == 0)
			{
				print("problem  as.treeshape")
				treebaselist$problematicwithas.treeshape <- c(treebaselist$problematicwithas.treeshape, i)
			}
			else
			{
				treebaselist$kind <- c(		treebaselist$kind , treebase[[i]]$kind)
				treebaselist$Nnode <- c(		treebaselist$Nnode, treebase[[i]]$Nnode )
				treebaselist$realleaves <- c(		treebaselist$realleaves, length(treebase[[i]]$tip.label))
				treebaselist$colless <- c(		treebaselist$colless, 	colless(as.treeshape(treebase[[i]])))	
				treebaselist$sackin <- c(		treebaselist$sackin, 	sackin(as.treeshape(treebase[[i]])))
				write.nexus(		treebase[[i]] , file=paste("treebaseforbeta", numbtrees , ".nex", sep=""))
				forbeta <- read.nexus(paste("treebaseforbeta", numbtrees , ".nex", sep=""))
				treebaselist$maxlik.betasplit <- c(		treebaselist$maxlik.betasplit, maxlik.betasplit(forbeta, up=10)$max_lik )				
			} 
			}
		}
		}	
	}
	else
	{
		{
		if ( treebase[[i]] == 0)
		{
			treebaselist$extinct <- c(treebaselist$extinct, 1)
		}
		else
		{
			treebaselist$samesp <- c(treebaselist$samesp, 1)
		}
		}
	}
	}
}

save(treebaselist, treeswithlength, file="treebase_output.RData")



#####Example on how to plot some trees
#par(mfrow=c(5,5))
#for (ptree in 1:25)
#{
#	plot(treebase[[ptree]])
#}

#par(mfrow=c(1,1))
#plot(treebase[[25]])

#drop.fossil(as.treeshape(treebase[[25]]))
#testeold <- new2old.phylo(treebase[[25]])
#class(testeold)
#testenew <- old2new.phylo(treebase[[25]])
#class(testenew)
#drop.fossil(testeold)
#prune.extinct.taxa(treebase[[25]])
#drop.fossil(treebase[[25]])
#class(treebase[[25]])