


# These two functions are where the hard-coded color and symbol options reside for PCA plotting

colorSelection = function(n, alpha = 255) {
  rawcolors = c('turquoise', 'red', 'yellow2', 'darkblue', 'limegreen', 'hotpink', 'blue', 
             'purple', 'brown', 'seagreen', 'darkorange', 'pink', 'firebrick4', 'olivedrab1', 'black')
  colorvector = col2rgb(rawcolors)
  colors = apply(colorvector, 2, function(x) rgb(x[1], x[2], x[3], alpha = alpha, maxColorValue = 255))
  
  if (n > length(colors)) { 
    warning("exceeding the maximum number of distinct colors; colors are recycled", immediate. = TRUE)
  }
  
  return(rep(colors, ceiling(n/length(colors)))[1:n])
}


pchSelection = function(n) {
  pch = c(15, 16, 17, 18, 1, 7, 8, 10, 2, 4, 12, 6, 3, 5)
  
  if (n > length(pch)) { 
    warning("exceeding the maximum number of distinct symbols; symbols are recycled", immediate. = TRUE)
  } 

  return(rep(pch, ceiling(n/length(pch)))[1:n])
}




# Script for calculating metrics on all simulated phylogenies and conducting PCA
# Function that conducts PCA on treeOutput set of specified tree metrics
# -- vars is a vector of column names to include, default is all vars
# -- models is a vector of model abbreviations

treeMetricsPCA = function(treeOutput, models = 'all', vars = 'all') {
  
  if (models[1] != 'all') {
    treeOutput = filter(treeOutput, model %in% models)
  }
  
  if (vars[1] == 'all') {
    vars = names(treeOutput[, 3:ncol(treeOutput)])
  }
  outputSubset = treeOutput[, c("model", "simID", vars)]
  outputSubsetNoNAs = na.omit(outputSubset)  
  
  pc = princomp(outputSubsetNoNAs[, names(outputSubsetNoNAs) %in% vars], cor = TRUE)
  pcaScores = cbind(outputSubsetNoNAs[, c("model", "simID")], pc$scores) 
  
  return(list(pcaScores = pcaScores, pcaLoadings = pc$loadings))
}


# Function for exploring tree shape as a function of between-model classifications

betweenModelPCAPlot = function(pcaOutput,          # dataframe with model, simID, and PC scores
                               xscore = 1,        # PC score plotted on the x-axis
                               yscore = 2,        # PC score plotted on the y-axis
                               colorBy = 'model', # any variable/parameter name by which to color points
                               pchBy = 'model',         # categorical variable/parameter name by which to specify point shape 
                               alpha = 255,       # transparency of symbols, where 255 is solid and 0 is totally transparent
                               ...) {
  
  # Reads in model classification codes and assigns colors and pchs
  require(gsheet)
  url = "https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=2047946073"
  
  modelClassification = gsheet2tbl(url)
  
  pcaScores = pcaOutput$pcaScores
  pcaLoadings = pcaOutput$pcaLoadings
  
  colorCode = data.frame(val = unique(modelClassification[, colorBy]), 
                         color = colorSelection(nrow(unique(modelClassification[, colorBy])), alpha))
  colorCode$color = as.character(colorCode$color)
  names(colorCode)[1] = colorBy
  
  pchCode = data.frame(val = unique(modelClassification[, pchBy]),
                       pch = pchSelection(nrow(unique(modelClassification[, pchBy]))))
  names(pchCode)[1] = pchBy
  
  plotOutput = left_join(pcaScores, modelClassification, by = c("model", "simID")) %>%
    left_join(colorCode, by = unname(colorBy)) %>%
    left_join(pchCode, by = unname(pchBy))

  maxAbsX = max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")])))
  maxAbsY = max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")])))
  
  plot(plotOutput[,paste("Comp.", xscore, sep = "")], plotOutput[, paste("Comp.", yscore, sep = "")], 
       col = plotOutput$color,  
       pch = plotOutput$pch, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-maxAbsY, maxAbsY),
       xlim = c(-maxAbsX, maxAbsX), ...)
  legend("topleft", 
         legend = c(toupper(colorBy), colorCode[,1]), bty = "n",
         col = c('white', colorCode[,2]), pch = 16, pt.cex = 2)
  
  legend("topright", 
         legend = c(toupper(pchBy), pchCode[,1]), bty = "n",
         col = c('white', rep('black', nrow(pchCode))), pch = c(16, pchCode[,2]), pt.cex = 2)
  # Only print the top variable loadings since it otherwise becomes messy
  mainLoadings = pcaLoadings[order(pcaLoadings[, xscore], decreasing = TRUE) <= 4 | 
                               order(pcaLoadings[,yscore], decreasing = TRUE) <= 4, 
                             c(xscore, yscore)]
  text(mainLoadings[, 1]*maxAbsX, mainLoadings[, 2]*maxAbsY, row.names(mainLoadings))  
}




# Function for exploring tree shape as a function of within-model parameters

withinModelPCAPlot = function(pcaOutput,          # dataframe with model, simID, and PC scores
                              modelAbbrev,             # specify the abbreviation of the model to explore 
                              xscore = 1,        # PC score plotted on the x-axis
                              yscore = 2,        # PC score plotted on the y-axis
                              colorBy,        # any variable/parameter name by which to color points
                              pchBy,            # categorical variable/parameter name by which to specify point shape 
                              alpha = 255,      # transparency of symbols where 255 is solid and 0 is totally transparent 
                              ...) {

  pcaScores = pcaOutput$pcaScores %>%
    filter(model == modelAbbrev)
  pcaLoadings = pcaOutput$pcaLoadings
  
  modelParams = read.csv(paste('parameters/', modelAbbrev, '_parameters.csv', sep = ''), header = TRUE)
  
  modelScores = left_join(pcaScores, modelParams, by = "simID")
  
  if (is.numeric(modelScores[, colorBy]) & length(unique(modelScores[, colorBy])) > 4) {
    
    shades <- rainbow(130)[100:1]
    percents <- as.integer(cut(modelScores[, colorBy], 100, include.lowest = TRUE, ordered = TRUE))
    modelScores$color = shades[percents]
    modelOutput = modelScores
    
  } else {
    
    colorCode = data.frame(val = unique(modelScores[, colorBy]), 
                           color = colorSelection(length(unique(modelScores[, colorBy])), alpha))
    colorCode$color = as.character(colorCode$color)
    names(colorCode)[1] = colorBy
    
    modelOutput = left_join(modelScores, colorCode, by = unname(colorBy))
  }
  
  if (class(pchBy) == 'numeric') {
    warning("Symbol size is best used for visualizing categorical variables", immediate. = TRUE)
  } else {
    pchCode = data.frame(val = unique(modelScores[, pchBy]),
                         pch = pchSelection(length(unique(modelScores[, pchBy]))))
    names(pchCode)[1] = pchBy
  }
  
  plotOutput = left_join(modelOutput, pchCode, by = unname(pchBy))
  
  maxAbsX = max(abs(range(plotOutput[,paste("Comp.", xscore, sep = "")])))
  maxAbsY = max(abs(range(plotOutput[,paste("Comp.", yscore, sep = "")])))
  
  plot(plotOutput[,paste("Comp.", xscore, sep = "")], plotOutput[, paste("Comp.", yscore, sep = "")], 
       col = plotOutput$color, 
       pch = plotOutput$pch, xlab = paste("PC", xscore), ylab = paste("PC", yscore), 
       ylim = c(-maxAbsY, maxAbsY),
       xlim = c(-maxAbsX, maxAbsX), ...)
  
  # color legend
  if (is.numeric(modelScores[, colorBy]) & length(unique(modelScores[, colorBy])) > 4) {
    maxvar = max(modelScores[, colorBy], na.rm = TRUE)
    minvar = min(modelScores[, colorBy], na.rm = TRUE)
    inc <- (maxvar - minvar) / 4
    legend.text <- round(c(minvar, minvar + inc, minvar + 2 * inc, minvar + 3 * inc, maxvar))
    
    legend("topleft", 
           legend = c(toupper(colorBy), legend.text), bty = "n",
           fill = c('white', shades[c(1, 25, 50, 75, 100)]))

  } else {

    legend("topleft", legend = c(toupper(colorBy), as.character(colorCode[,1])), bty = "n",
           col = c('white', as.character(colorCode[,2])), pch = 16, pt.cex = 2)
  }

  legend("topright", 
         legend = c(toupper(pchBy), as.character(pchCode[,1])), bty = "n",
         col = c('white', rep('black', nrow(pchCode))), pch = c(16, pchCode[,2]), pt.cex = 2)
  
  # Only print the top variable loadings since it otherwise becomes messy
  mainLoadings = pcaLoadings[order(pcaLoadings[, xscore], decreasing = TRUE) <= 4 | 
                               order(pcaLoadings[,yscore], decreasing = TRUE) <= 4, 
                             c(xscore, yscore)]
  text(mainLoadings[, 1]*maxAbsX, mainLoadings[, 2]*maxAbsY, row.names(mainLoadings))  
}




# Function for exploring tree shape as a function of between-model classifications, plotting raw
# tree metrics rather than PC scores

betweenModelVarPlot = function(treeOutput,          # dataframe with model, simID, and PC scores
                               xvar = 'Beta',        # tree metric to be plotted on the x-axis
                               yvar = 'Gamma',        # tree metric to be plotted on the y-axis
                               colorBy = 'model', # any variable/parameter name by which to color points
                               pchBy = 'model',         # categorical variable/parameter name by which to specify point shape 
                               alpha = 255,      # transparency of symbols where 255 is solid and 0 is totally transparent 
                               ...) {
  
  # Reads in model classification codes and assigns colors and pchs
  require(gsheet)
  url = "https://docs.google.com/spreadsheets/d/1pcUuINauW11cE5OpHVQf_ZuzHzhm2VJkCn7-lSEJXYI/edit#gid=2047946073"
  
  modelClassification = gsheet2tbl(url)
  
  colorCode = data.frame(val = unique(modelClassification[, colorBy]), 
                         color = colorSelection(nrow(unique(modelClassification[, colorBy])), alpha))
  colorCode$color = as.character(colorCode$color)
  names(colorCode)[1] = colorBy
  
  pchCode = data.frame(val = unique(modelClassification[, pchBy]),
                       pch = pchSelection(nrow(unique(modelClassification[, pchBy]))))
  names(pchCode)[1] = pchBy
  
  plotOutput = left_join(treeOutput, modelClassification, by = c("model", "simID")) %>%
    left_join(colorCode, by = unname(colorBy)) %>%
    left_join(pchCode, by = unname(pchBy))
  
  miny = min(plotOutput[, yvar], na.rm = TRUE)
  maxy = max(plotOutput[, yvar], na.rm = TRUE)
  minx = min(plotOutput[, xvar], na.rm = TRUE)
  maxx = max(plotOutput[, xvar], na.rm = TRUE)
  
  plot(plotOutput[, xvar], plotOutput[, yvar], 
       col = plotOutput$color,  
       pch = plotOutput$pch, xlab = xvar, ylab = yvar, 
       xlim = c(minx - (maxx - minx)/8, maxx + (maxx - minx)/8), 
       ...)
  legend("topleft", 
         legend = c(toupper(colorBy), colorCode[,1]), bty = "n",
         col = c('white', colorCode[,2]), pch = 16, pt.cex = 2)
  
  legend("topright", 
         legend = c(toupper(pchBy), pchCode[,1]), bty = "n",
         col = c('white', rep('black', nrow(pchCode))), pch = c(16, pchCode[,2]), pt.cex = 2)
}


# Function for exploring tree shape metrics as a function of within-model parameters

withinModelVarPlot = function(treeOutput,          # dataframe with model, simID, and tree shape metrics
                              modelAbbrev,             # specify the abbreviation of the model to explore 
                              modelParams = NULL,      # dataframe with simID and columns for all model parameters,
                                                  #   if NULL, then read it in from 'parameters' folder
                              xvar = 'Beta',        # PC score plotted on the x-axis
                              yvar = 'Gamma',        # PC score plotted on the y-axis
                              colorBy,        # any variable/parameter name by which to color points
                              pchBy,            # categorical variable/parameter name by which to specify point shape 
                              alpha = 255,      # transparency of symbols where 255 is solid and 0 is totally transparent 
                              ...) {
  
  if (is.null(modelParams)) {
    modelParams = read.csv(paste('parameters/', modelAbbrev, '_parameters.csv', sep = ''), header = TRUE, stringsAsFactors = FALSE)
  }
  modelOutput = filter(treeOutput, model == modelAbbrev) %>%
    left_join(modelParams, by = c("model", "simID"))
  
  if (is.numeric(modelOutput[, colorBy]) & length(unique(modelOutput[, colorBy])) > 4) {
    
    shades <- rainbow(130)[100:1]
    percents <- as.integer(cut(modelOutput[, colorBy], 100, include.lowest = TRUE, ordered = TRUE))
    modelOutput$color = shades[percents]
   
  } else {
    
    colorCode = data.frame(val = unique(modelOutput[, colorBy]), 
                           color = colorSelection(length(unique(modelOutput[, colorBy])), alpha))
    colorCode$color = as.character(colorCode$color)
    if (class(colorCode$val) == 'factor') { colorCode$val = as.character(colorCode$val) }
    names(colorCode)[1] = colorBy
    
    modelOutput = left_join(modelOutput, colorCode, by = unname(colorBy))
  }
  
  if (class(pchBy) == 'numeric' & length(unique(modelOutput[, pchBy])) > 5) {
    warning("Symbol size is best used for visualizing categorical variables", immediate. = TRUE)
  } else {
    pchCode = data.frame(val = unique(modelOutput[, pchBy]),
                         pch = pchSelection(length(unique(modelOutput[, pchBy]))))
    if (class(pchCode$val) == 'factor') { pchCode$val = as.character(pchCode$val) }
    names(pchCode)[1] = pchBy
  }
  
  plotOutput = left_join(modelOutput, pchCode, by = unname(pchBy))
  
  miny = min(plotOutput[, yvar], na.rm = TRUE)
  maxy = max(plotOutput[, yvar], na.rm = TRUE)
  minx = min(plotOutput[, xvar], na.rm = TRUE)
  maxx = max(plotOutput[, xvar], na.rm = TRUE)
  
  plot(plotOutput[, xvar], plotOutput[, yvar], 
       col = plotOutput$color, 
       pch = plotOutput$pch, xlab = xvar, ylab = yvar, 
       xlim = c(minx - (maxx - minx)/8, maxx + (maxx - minx)/8), 
       ...)
  
  # color legend
  if (is.numeric(modelOutput[, colorBy]) & length(unique(modelOutput[, colorBy])) > 4) {
    maxvar = max(modelOutput[, colorBy], na.rm = TRUE)
    minvar = min(modelOutput[, colorBy], na.rm = TRUE)
    inc <- (maxvar - minvar) / 4
    legend.text <- round(c(minvar, minvar + inc, minvar + 2 * inc, minvar + 3 * inc, maxvar))
    
    legend("topleft", 
           legend = c(toupper(colorBy), legend.text), bty = "n",
           fill = c('white', shades[c(1, 25, 50, 75, 100)]))
    
  } else {
    
    legend("topleft", legend = c(toupper(colorBy), as.character(colorCode[,1])), bty = "n",
           col = c('white', as.character(colorCode[,2])), pch = 16, pt.cex = 2)
  }
  
  legend("topright", 
         legend = c(toupper(pchBy), as.character(pchCode[,1])), bty = "n",
         col = c('white', rep('black', nrow(pchCode))), pch = c(16, pchCode[,2]), pt.cex = 2)
  
}
