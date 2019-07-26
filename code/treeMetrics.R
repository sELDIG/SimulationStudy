# Functions for calculating various metrics on phylogenies




# This function is a modification of the maxlik.betasplit() function
# in the apTreeshape package for calculating 'beta'

maxlik.betasplit.AH = function (phylo, up = 10, remove.outgroup = FALSE, confidence.interval = "none", 
                                conf.level = 0.95, size.bootstrap = 100) 
{
  vrais.aldous.fin <- function(i, n, b) {
    # Code commented out below is from the original maxlik.betasplit function.
    # Due to underflow errors, beta was switched out for lbeta in the 
    # manner of Purvis et al. 2011, Phil Trans Roy Soc B 366: 2462-2477
    # (Purvis, pers. comm.)
    #aux <- beta(b + i + 1, b + n - i + 1)/beta(i + 1, n - 
    #                                             i + 1)
    #if (is.na(aux) | (aux == Inf)) 
    #  aux <- (i/n)^b * (1 - i/n)^(b)
    aux <- exp(lbeta(b + i + 1, b + n - i + 1) - lbeta(i + 1, n - i + 1))
    return(aux)
  }
  bbalance <- function(phylo) {
    return(t(apply(balance(phylo), FUN = function(x) {
      c(x[1], x[1] + x[2])
    }, MARGIN = 1)))
  }
  renorm.aldous <- function(n, beta) {
    return(sum(sapply(1:(n - 1), FUN = vrais.aldous.fin, 
                      n = n, b = beta)))
  }
  vrais.aldous.renorm <- function(i, n, beta) {
    return(vrais.aldous.fin(i, n, beta)/renorm.aldous(n, 
                                                      beta))
  }
  logvrais.aldous.phylo <- function(b, phylo, remove.outgroup = TRUE) {
    if (class(phylo) == "treeshape") 
      bal <- bbalance(as.phylo(phylo))
    if (class(phylo) == "phylo") 
      bal <- bbalance(phylo)
    if (remove.outgroup) {
      if ((bal[1, 1] <= 2) || ((bal[1, 2] - bal[1, 1]) <= 
                               2)) 
        bal <- bal[-1, ]
    }
    return(sum(log(apply(bal, FUN = function(x, b) {
      return(vrais.aldous.renorm(x[1], x[2], b))
    }, b = b, MARGIN = 1))))
  }
  logvrais.aldous.bal <- function(b, bal) {
    return(sum(log(apply(bal, FUN = function(x, b) {
      return(vrais.aldous.renorm(x[1], x[2], b))
    }, b = b, MARGIN = 1))))
  }
  if (class(phylo) == "treeshape") 
    bal <- bbalance(as.phylo(phylo))
  if (class(phylo) == "phylo") 
    bal <- bbalance(phylo)
  if (class(phylo) != "phylo" && class(phylo) != "treeshape") {
    print("The phylogeny shall be of class phylo or treeshape")
    return
  }
  if (remove.outgroup) {
    if ((bal[1, 1] <= 2) || ((bal[1, 2] - bal[1, 1]) <= 2)) 
      bal <- bal[-1, ]
  }
  nb.tip <- max(bal[1, ])
  optim_lik_aldous <- function(phylo, remove.outgroup) {
    optimize(f = function(x) {
      logvrais.aldous.phylo(x, phylo, remove.outgroup)
    }, lower = -2, upper = up, maximum = TRUE)
  }
  optim_lik_aldous_bal <- function(bal) {
    optimize(f = function(x) {
      logvrais.aldous.bal(x, bal)
    }, lower = -2, upper = up, maximum = TRUE)
  }
  res <- optim_lik_aldous(phylo, remove.outgroup)
  if (confidence.interval == "bootstrap") {
    nb_b <- dim(bal)[1]
    fun_aux <- function() {
      optim_lik_aldous_bal(bal[sample(1:nb_b, size = nb_b, 
                                      replace = TRUE), ])$maximum
    }
    thebeta <- replicate(size.bootstrap, fun_aux())
    up.conf <- 1 - ((1 - conf.level)/2)
    low.conf <- (1 - conf.level)/2
    conf_interval <- quantile(thebeta, c(low.conf, up.conf))
  }
  if (confidence.interval == "profile") {
    if ((res$objective - 1.92) - logvrais.aldous.phylo(-2, 
                                                       phylo, remove.outgroup) < 0) 
      low.conf <- (-2)
    else low.conf <- uniroot(f = function(x) {
      (res$objective - 1.92) - logvrais.aldous.phylo(x, 
                                                     phylo, remove.outgroup)
    }, lower = -2, upper = res$maximum)$root
    if ((res$objective - 1.92) - logvrais.aldous.phylo(up, 
                                                       phylo, remove.outgroup) < 0) 
      up.conf <- up
    else up.conf <- uniroot(f = function(x) {
      (res$objective - 1.92) - logvrais.aldous.phylo(x, 
                                                     phylo, remove.outgroup)
    }, lower = res$maximum, upper = up)$root
    conf_interval <- c(low.conf, up.conf)
  }
  if (confidence.interval == "none") {
    conf_interval <- NULL
  }
  return(list(max_lik = res$maximum, conf_interval = conf_interval))
}


# Function for measuring various attributes on a phylogeny in Newick format

treeMetrics = function(treeInput) {

  if(class(treeInput) != "phylo") {
    stop("The object passed to treeMetrics() is not of class 'phylo'")
    }
  
  require(ape)
  require(caper)
  require(geiger)
  require(picante)
  require(apTreeshape)
  require(dispRity)
  require(RPANDA)
  require(nLTT)
  require(stringr)
  require(Rfast)
  
  # Drop root edge
  treeInput$root.edge = 0

  # Prune out extinct species
  tree = tryCatch({
    drop.extinct(treeInput)
  }, error = function(e) {
    tree = drop.fossil(treeInput)
  })
  
  # Richness
  S = length(tree$tip.label)
  
  # Absolute tree length, which is the diagonal of the vcv matrix
  v.matrix = vcv(tree, corr=F)
  tree.length = diag(v.matrix)[1]
  
  # Tree with branch lengths scaled by total tree length
  tree.scaled = tree
  tree.scaled$edge.length = tree$edge.length/tree.length
  
  # PD (phylogenetic diversity, summed branch lengths) on scaled tree
  PD = sum(tree.scaled$edge.length)
  
  # Create treeshape objects for certain calculations
  tree.shape = as.treeshape(tree)
  tree.shape.scaled = as.treeshape(tree.scaled)
  
  #Pybus & Harvey (2000)'s gamma statistic
  gamma.stat = gammaStat(tree.scaled)

  # Some metrics are difficult to calculate on very large trees.
  # Calculate Blum & Francois (2006)'s Beta metric of tree imbalance using apTreeshape package
  # Also calculate RPANDA spectral density metrics (Lewitus & Morlon 2016)
  
  # both analyses seem to bonk on very large phylogenies, so only try calculating for fewer than 6000 species
  if(S < 6000) {
    # beta
    beta.out = maxlik.betasplit.AH(tree.scaled)
    beta.stat = beta.out$max_lik
    
    # RPANDA
    MGL = spectR(tree)
    MGL_principal_eigenvalue = MGL$principal_eigenvalue 
    MGL_asymmetry = MGL$asymmetry  
    MGL_peakedness = MGL$peakedness
    MGL_eigengap = MGL$eigengap
    
    # Mean Pairwise Distance (in scaled tree)
    pairwise.dists = dist.nodes(tree.scaled)
    sum.pairwise.dists = upper_tri(pairwise.dists, diag = FALSE, suma = TRUE)
    MPD = sum.pairwise.dists/(ncol(pairwise.dists)*(ncol(pairwise.dists)-1)/2)
    
    # Variance Pairwise Distance (in scaled tree)
    VPD = tryCatch({
      var(as.vector(tri.pairwise.dists))
    }, error = function(e) {
      VPD = NA
    })
    

  } else {
    beta.stat = NA
    MGL_principal_eigenvalue = NA
    MGL_asymmetry = NA  
    MGL_peakedness = NA
    MGL_eigengap = NA
    MPD = NA
    VPD = NA
  }
  
  # Colless index
  Colless = colless(tree.shape)
  
  # Sackin index
  Sackin = sackin(tree.shape)
  
  # logarithm of the ratio of the likelihoods under the Yule model and the PDA model
  Yule.PDA.ratio = shape.statistic(tree.shape)
  
  # Mean Root Distance (MRD, also abbreviated N-bar, Shao & Sokal 1990) (Agapow & Purvis 2002)
  # Modified from code originally by ELIOT MILLER 25 AUGUST 2011
  phylo.bl1 <- compute.brlen(tree, 1)
  all.dist <- dist.nodes(phylo.bl1)
  root.dist <- all.dist[length(tree$tip.label)+1, 1:length(tree$tip.label)]
  MRD = mean(root.dist)
  
  # Variance in root distances across the tips, Shao & Sokal 1990 (Agapow & Purvis 2002)
  VRD = var(root.dist)
  
  # Calculate Helmus et al. 2007's Phylogenetic Species Variability (PSV) metric
  # based on the phylogenetic variance/covariance matrix
  # While Helmus et al.'s equation is as follows
  # PSV = (n*trace(C) - sum(C))/(n*(n-1)) 
  # this only holds for a tree with contemporaneous tips where C is the correlation matrix.
  # When the tree is not ultra-metric, Algar et al. 2009 (Ecol Lett) show the calculation
  # using the variance-covariance matrix, V:
  # PSV = (n*trace(V) - sum(V))/(trace(V)*(n-1)) 
  
  n = nrow(v.matrix)
  PSV = (n*sum(diag(v.matrix)) - sum(v.matrix))/(sum(diag(v.matrix))*(n-1))
  
  # mean I' from Purvis et al. 2002. Evaluating phylogenetic tree shape: two modifications to Fusco & Cronk's method
  if(tree$Nnode >= 3) {
    tree$node.label = NULL
    fusco = caper::fusco.test(tree)
    mean.Iprime = fusco$mean.Iprime
  } else {
    mean.Iprime = NA
  }
  
  # nLTT statistic on first example tree
  utils::data(exampleTrees, package = "nLTT") # Set of birth-death trees
  nLTT_stat <- nLTT::nLTTstat(
    tree1 = exampleTrees[[1]], # Can be changed to generated Yule tree or any other tree
    tree2 = treeInput
  )

  
  return(list(S = S, tree.length = tree.length, PD = PD, gamma = gamma.stat, beta = beta.stat, 
              Colless = Colless, Sackin = Sackin, Yule.PDA.ratio = Yule.PDA.ratio, MRD = MRD, 
              VRD = VRD, PSV = PSV, mean.Iprime = mean.Iprime,
              MPD = MPD, VPD = VPD, 
              MGL_principal_eigenvalue = MGL_principal_eigenvalue, MGL_asymmetry = MGL_asymmetry, 
              MGL_peakedness = MGL_peakedness, MGL_eigengap = MGL_eigengap, nLTT_stat = nLTT_stat))
}


# Run treeMetrics for many trees and save output as it goes

metricsForManyTrees = function(treefiles = NULL, minimumTreeSize = 20, fileOut, append = TRUE) {
  
  if(is.null(treefiles)) {
    treefiles = list.files('trees')[grepl(".tre", list.files("trees"))]
  }  
  
  if(!append) {
    
    treeOutput = data.frame(model = NA, simID = NA, S = NA, tree.length = NA, PD = NA, gamma = NA, 
                            beta = NA, Colless = NA, Sackin = NA, Yule.PDA.ratio = NA, MRD = NA, 
                            VRD = NA, PSV = NA, mean.Iprime = NA, MPD = NA, VPD = NA, 
                            MGL_principal_eigenvalue = NA, MGL_asymmetry = NA,
                            MGL_peakedness = NA, MGL_eigengap = NA, nLTT_stat = NA)
    
    write.table(treeOutput, fileOut, sep = '\t', row.names = FALSE)
    
  }

  for (treefile in treefiles) {
    
    tree = read.tree(paste("trees/", treefile, sep = ""))
    
    if(tree$Nnode + 1 >= minimumTreeSize) {
      model = str_extract(treefile, "^[A-Za-z]*")
      simID = str_extract(treefile, "[0-9]+")
      
      print(paste(treefile, Sys.time()))
      
      metrics = tryCatch({
        treeMetrics(tree)
      }, error = function(e) {
         metrics =  data.frame(model = NA, simID = NA, S = NA, tree.length = NA, PD = NA, gamma = NA, 
                               beta = NA, Colless = NA, Sackin = NA, Yule.PDA.ratio = NA, MRD = NA, 
                               VRD = NA, PSV = NA, mean.Iprime = NA, MPD = NA, VPD = NA, 
                               MGL_principal_eigenvalue = NA, MGL_asymmetry = NA,
                               MGL_peakedness = NA, MGL_eigengap = NA, nLTT_stat = NA)
      })
      
      sink(fileOut, append = TRUE)
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
      

    } else {
      print(paste(treefile, "skipped -- not enough species"))
    }
    
  }  
}

  
