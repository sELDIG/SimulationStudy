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
  
  # Prune out extinct species
  tree = drop.extinct(treeInput)
  tree.shape = as.treeshape(tree)
  
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
  
  # Richness
  S = length(tree$tip.label)
    
  #Pybus & Harvey (2000)'s gamma statistic
  gamma.stat = gammaStat(tree)
  
  #Calculate Blum & Francois (2006)'s Beta metric of tree imbalance using apTreeshape package
  # --seems to bonk on very large phylogenies, so only try calculating for fewer than 6000 species
  if(S < 6000) {
    beta.out = maxlik.betasplit.AH(tree)
    beta.stat = beta.out$max_lik
  } else {
    beta.stat = NA
  }
  
  # Colless index
  Colless = colless(tree.shape)
  
  # Sackin index
  Sackin = sackin(tree.shape)
  
  # logarithm of the ratio of the likelihoods under the Yule model and the PDA model
  shape.stat = shape.statistic(tree.shape)
  
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
  
  v.matrix = vcv(tree, corr=F)
  n = nrow(v.matrix)
  PSV = (n*sum(diag(v.matrix)) - sum(v.matrix))/(sum(diag(v.matrix))*(n-1))
  
  # mean I' from Purvis et al. 2002. Evaluating phylogenetic tree shape: two modifications to Fusco & Cronk's method
  if(tree$Nnode >= 3) {
    fusco = caper::fusco.test(tree)
    mean.Iprime = fusco$mean.Iprime
  } else {
    mean.Iprime = NA
  }
  
  # TO DO: 
  # --identify additional metrics
  # --incoporate RPANDA, ClaDS, etc. output.
  
  
  
  
  return(list(S = S, gamma = gamma.stat, beta = beta.stat, Colless = Colless, 
              Sackin = Sackin, shape = shape.stat, MRD = MRD, VRD = VRD, PSV = PSV, mean.Iprime = mean.Iprime))
}


  
  
